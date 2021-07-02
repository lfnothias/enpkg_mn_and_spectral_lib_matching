import pandas as pd
import numpy as np
import zipfile
import glob
import os
import sys
import time
import shlex
import subprocess
import contextlib
import io
from tqdm import tqdm
from tqdm import tqdm_notebook
from tqdm.contrib import tzip
from opentree import OT
import json
from pandas import json_normalize
import yaml
from pandas import json_normalize
import networkx as nx

from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import select_by_intensity
from matchms.filtering import select_by_mz
from matchms.similarity import PrecursorMzMatch
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from matchms.similarity import ModifiedCosine

# We deactivate the iloc warning see https://stackoverflow.com/a/20627316
pd.options.mode.chained_assignment = None  # default='warn'

with open (r'configs/default/default.yaml') as file:    
    params_list = yaml.load(file, Loader=yaml.FullLoader)

repository_path = params_list['paths'][0]['repository_path']
spectra_suffix = params_list['paths'][1]['spectra_suffix']
metadata_sample_suffix = params_list['paths'][2]['metadata_sample_suffix']
metadata_path = params_list['paths'][3]['metadata_path']
db_file_path = params_list['paths'][4]['db_file_path']
adducts_pos_path = params_list['paths'][5]['adducts_pos_path']
adducts_neg_path = params_list['paths'][6]['adducts_neg_path']

parent_mz_tol = params_list['spectral_match_params'][0]['parent_mz_tol']
msms_mz_tol = params_list['spectral_match_params'][1]['msms_mz_tol']
min_cos = params_list['spectral_match_params'][2]['min_cos']
min_peaks = params_list['spectral_match_params'][3]['min_peaks']
match_score = params_list['spectral_match_params'][4]['match_score']

mn_parent_mz_tol = params_list['networking_params'][0]['mn_parent_mz_tol']
mn_msms_mz_tol = params_list['networking_params'][1]['mn_msms_mz_tol']
mn_score_cutoff = params_list['networking_params'][2]['mn_score_cutoff']
mn_max_links = params_list['networking_params'][3]['mn_max_links']
mn_top_n = params_list['networking_params'][4]['mn_top_n']
mn_score = params_list['networking_params'][5]['mn_score']

top_to_output= params_list['repond_params'][0]['top_to_output']
ppm_tol = params_list['repond_params'][1]['ppm_tol']
polarity = params_list['repond_params'][2]['polarity']
organism_header = params_list['repond_params'][3]['organism_header']
sampletype_header = params_list['repond_params'][4]['sampletype_header']
use_post_taxo = params_list['repond_params'][5]['use_post_taxo']
top_N_chemical_consistency = params_list['repond_params'][6]['top_N_chemical_consistency']
min_score_ms1 = params_list['repond_params'][7]['min_score_ms1']

isdb_output_suffix = params_list['output_params'][0]['isdb_output_suffix']
mn_output_suffix = params_list['output_params'][1]['mn_output_suffix']
repond_table_suffix = params_list['output_params'][2]['repond_table_suffix']

# Defining functions 

@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = io.StringIO()
    yield
    sys.stdout = save_stdout

def metadata_processing(spectrum):
    spectrum = default_filters(spectrum)
    return spectrum


def peak_processing(spectrum):
    spectrum = default_filters(spectrum)
    spectrum = normalize_intensities(spectrum)
    spectrum = select_by_intensity(spectrum, intensity_from=0.01)
    spectrum = select_by_mz(spectrum, mz_from=10, mz_to=1000)
    return spectrum

def connected_component_subgraphs(G):
            for c in nx.connected_components(G):
                yield G.subgraph(c)


######################## MN FUNCTIONS FROM MATCHMS ####################################3

from typing import Tuple
from matchms import Scores               
def get_top_hits(scores: Scores, identifier_key: str = "spectrumid",
                 top_n: int = 25, search_by: str = "queries",
                 ignore_diagonal: bool = False) -> Tuple[dict, dict]:
    """Get top_n highest scores (and indices) for every entry.
    Parameters
    ----------
    scores
        Matchms Scores object containing all similarities.
    identifier_key
        Metadata key for unique intentifier for each spectrum in scores.
        Will also be used for the naming the network nodes. Default is 'spectrumid'.
    top_n
        Return the indexes and scores for the top_n highest scores. Scores between
        a spectrum with itself (diagonal of scores.scores) will not be taken into
        account.
    search_by
        Chose between 'queries' or 'references' which decides if the top_n matches
        for every spectrum in scores.queries or in scores.references will be
        collected and returned
    ignore_diagonal
        Set to True if scores.scores is symmetric (i.e. if references and queries
        were the same) and if scores between spectra with themselves should be
        excluded.
    """
    assert search_by in ["queries", "references"], \
        "search_by must be 'queries' or 'references"

    similars_idx = dict()
    similars_scores = dict()

    if search_by == "queries":
        for i, spec in enumerate(scores.queries):
            spec_id = spec.get(identifier_key)
            idx = scores.similarity_function.sort(scores.scores[:, i])
            if ignore_diagonal:
                similars_idx[spec_id] = idx[idx != i][:top_n]
            else:
                similars_idx[spec_id] = idx[:top_n]
            similars_scores[spec_id] = scores.scores[similars_idx[spec_id], i]
    elif search_by == "references":
        for i, spec in enumerate(scores.references):
            spec_id = spec.get(identifier_key)
            idx = scores.similarity_function.sort(scores.scores[i, :])
            if ignore_diagonal:
                similars_idx[spec_id] = idx[idx != i][:top_n]
            else:
                similars_idx[spec_id] = idx[:top_n]
            similars_scores[spec_id] = scores.scores[i, similars_idx[spec_id]]
    return similars_idx, similars_scores

def create_network(self, scores: Scores):
    """
    Function to create network from given top-n similarity values. Expects that
    similarities given in scores are from an all-vs-all comparison including all
    possible pairs.
    Parameters
    ----------
    scores
        Matchms Scores object containing all spectrums and pair similarities for
        generating a network.
    """
    assert self.top_n >= self.max_links, "top_n must be >= max_links"
    assert np.all(scores.queries == scores.references), \
        "Expected symmetric scores object with queries==references"
    unique_ids = list({s.get(self.identifier_key) for s in scores.queries})

    # Initialize network graph, add nodes
    msnet = nx.Graph()
    msnet.add_nodes_from(unique_ids)

    # Collect location and score of highest scoring candidates for queries and references
    similars_idx, similars_scores = get_top_hits(scores, identifier_key=self.identifier_key,
                                                    top_n=self.top_n,
                                                    search_by="queries",
                                                    ignore_diagonal=True)
    similars_scores = self._select_edge_score(similars_scores, scores.scores.dtype)

    # Add edges based on global threshold (cutoff) for weights
    for i, spec in enumerate(scores.queries):
        query_id = spec.get(self.identifier_key)

        ref_candidates = np.array([scores.references[x].get(self.identifier_key)
                                        for x in similars_idx[query_id]])
        idx = np.where((similars_scores[query_id] >= self.score_cutoff) &
                            (ref_candidates != query_id))[0][:self.max_links]
        if self.link_method == "single":
            new_edges = [(query_id, str(ref_candidates[x]),
                            float(similars_scores[query_id][x])) for x in idx]
        elif self.link_method == "mutual":
            new_edges = [(query_id, str(ref_candidates[x]),
                            float(similars_scores[query_id][x]))
                            for x in idx if i in similars_idx[ref_candidates[x]][:]]
        else:
            raise ValueError("Link method not kown")

        msnet.add_weighted_edges_from(new_edges)

    if not self.keep_unconnected_nodes:
        msnet.remove_nodes_from(list(nx.isolates(msnet)))
    self.graph = msnet


from typing import Optional
import networkx as nx
import numpy
from matchms import Scores


class SimilarityNetwork:
    """Create a spectal network from spectrum similarities.
    For example
    .. testcode::
        import numpy as np
        from matchms import Spectrum, calculate_scores
        from matchms.similarity import ModifiedCosine
        from matchms.networking import SimilarityNetwork
        spectrum_1 = Spectrum(mz=np.array([100, 150, 200.]),
                              intensities=np.array([0.7, 0.2, 0.1]),
                              metadata={"precursor_mz": 100.0,
                                        "testID": "one"})
        spectrum_2 = Spectrum(mz=np.array([104.9, 140, 190.]),
                              intensities=np.array([0.4, 0.2, 0.1]),
                              metadata={"precursor_mz": 105.0,
                                        "testID": "two"})
        # Use factory to construct a similarity function
        modified_cosine = ModifiedCosine(tolerance=0.2)
        spectrums = [spectrum_1, spectrum_2]
        scores = calculate_scores(spectrums, spectrums, modified_cosine)
        ms_network = SimilarityNetwork(identifier_key="testID")
        ms_network.create_network(scores)
        nodes = list(ms_network.graph.nodes())
        nodes.sort()
        print(nodes)
    Should output
    .. testoutput::
        ['one', 'two']
    """
    def __init__(self, identifier_key: str = "spectrumid",
                 top_n: int = 20,
                 max_links: int = 10,
                 score_cutoff: float = 0.7,
                 link_method: str = 'single',
                 keep_unconnected_nodes: bool = True):
        """
        Parameters
        ----------
        identifier_key
            Metadata key for unique intentifier for each spectrum in scores.
            Will also be used for the naming the network nodes. Default is 'spectrumid'.
        top_n
            Consider edge between spectrumA and spectrumB if score falls into
            top_n for spectrumA or spectrumB (link_method="single"), or into
            top_n for spectrumA and spectrumB (link_method="mutual"). From those
            potential links, only max_links will be kept, so top_n must be >= max_links.
        max_links
            Maximum number of links to add per node. Default = 10.
            Due to incoming links, total number of links per node can be higher.
            The links are populated by looping over the query spectrums.
            Important side note: The max_links restriction is strict which means that
            if scores around max_links are equal still only max_links will be added
            which can results in some random variations (sorting spectra with equal
            scores restuls in a random order of such elements).
        score_cutoff
            Threshold for given similarities. Edges/Links will only be made for
            similarities > score_cutoff. Default = 0.7.
        link_method
            Chose between 'single' and 'mutual'. 'single will add all links based
            on individual nodes. 'mutual' will only add links if that link appears
            in the given top-n list for both nodes.
        keep_unconnected_nodes
            If set to True (default) all spectra will be included as nodes even
            if they have no connections/edges of other spectra. If set to False
            all nodes without connections will be removed.
        """
        # pylint: disable=too-many-arguments
        self.identifier_key = identifier_key
        self.top_n = top_n
        self.max_links = max_links
        self.score_cutoff = score_cutoff
        self.link_method = link_method
        self.keep_unconnected_nodes = keep_unconnected_nodes
        self.graph: Optional[nx.Graph] = None
        """NetworkX graph. Set after calling create_network()"""

    @staticmethod
    def _select_edge_score(similars_scores: dict, scores_type: numpy.dtype):
        """Chose one value if score contains multiple values (e.g. "score" and "matches")"""
        if len(scores_type) > 1 and "score" in scores_type.names:
            return {key: value["score"] for key, value in similars_scores.items()}
        if len(scores_type) > 1:  # Assume that first entry is desired score
            return {key: value[0] for key, value in similars_scores.items()}
        return similars_scores

    def create_network(self, scores: Scores):
        """
        Function to create network from given top-n similarity values. Expects that
        similarities given in scores are from an all-vs-all comparison including all
        possible pairs.
        Parameters
        ----------
        scores
            Matchms Scores object containing all spectrums and pair similarities for
            generating a network.
        """
        assert self.top_n >= self.max_links, "top_n must be >= max_links"
        assert numpy.all(scores.queries == scores.references), \
            "Expected symmetric scores object with queries==references"
        unique_ids = list({s.get(self.identifier_key) for s in scores.queries})

        # Initialize network graph, add nodes
        msnet = nx.Graph()
        msnet.add_nodes_from(unique_ids)

        # Collect location and score of highest scoring candidates for queries and references
        similars_idx, similars_scores = get_top_hits(scores, identifier_key=self.identifier_key,
                                                     top_n=self.top_n,
                                                     search_by="queries",
                                                     ignore_diagonal=True)
        similars_scores = self._select_edge_score(similars_scores, scores.scores.dtype)

        # Add edges based on global threshold (cutoff) for weights
        for i, spec in enumerate(scores.queries):
            query_id = spec.get(self.identifier_key)

            ref_candidates = numpy.array([scores.references[x].get(self.identifier_key)
                                          for x in similars_idx[query_id]])
            idx = numpy.where((similars_scores[query_id] >= self.score_cutoff) &
                              (ref_candidates != query_id))[0][:self.max_links]
            if self.link_method == "single":
                new_edges = [(query_id, str(ref_candidates[x]),
                              float(similars_scores[query_id][x])) for x in idx]
            elif self.link_method == "mutual":
                new_edges = [(query_id, str(ref_candidates[x]),
                              float(similars_scores[query_id][x]))
                             for x in idx if i in similars_idx[ref_candidates[x]][:]]
            else:
                raise ValueError("Link method not kown")

            msnet.add_weighted_edges_from(new_edges)

        if not self.keep_unconnected_nodes:
            msnet.remove_nodes_from(list(nx.isolates(msnet)))
        self.graph = msnet

    def export_to_graphml(self, filename: str):
        """Save the network as .graphml file.
        Parameters
        ----------
        filename
            Specify filename for exporting the graph.
        """
        if not self.graph:
            raise ValueError("No network found. Make sure to first run .create_network() step")
        nx.write_graphml_lxml(self.graph, filename)



















# timer is started
start_time = time.time()

print('Cleaning the spectral database metadata fields ...')

spectrums_db = list(load_from_mgf(db_file_path))


spectrums_db_cleaned = [metadata_processing(s) for s in spectrums_db]

print('A total of %s clean spectra were found in the spectral library.' % len(spectrums_db_cleaned))


if polarity == 'pos':
    adducts_df = pd.read_csv(adducts_pos_path, compression='gzip', sep='\t')
else:
    adducts_df = pd.read_csv(adducts_neg_path, compression='gzip', sep='\t')

adducts_df['min'] = adducts_df['adduct_mass'] - \
    int(ppm_tol) * (adducts_df['adduct_mass'] / 1000000)
adducts_df['max'] = adducts_df['adduct_mass'] + \
    int(ppm_tol) * (adducts_df['adduct_mass'] / 1000000)


repository_path_list = [x[0] for x in os.walk(repository_path)]
repository_path_list.remove(repository_path)

for sample_dir in repository_path_list:

    sample = sample_dir.split(repository_path,1)[1]
    #os.chdir(sample_dir)

    if len(glob.glob(sample_dir+'/*'+spectra_suffix)) != 0 :
        spectra_file_path = glob.glob(sample_dir+'/*'+ spectra_suffix)[0]

    if len(glob.glob(sample_dir+'/*'+metadata_sample_suffix)) != 0 :
        metadata_sample_path = glob.glob(sample_dir +'/*' + metadata_sample_suffix)[0]
        samples_metadata = pd.read_csv(metadata_sample_path, sep='\t')
    else:
        continue
    
    if samples_metadata[sampletype_header][0] == 'sample':
        print('Treating file ' + sample)

        isdb_results_path = sample_dir + '/' + sample + isdb_output_suffix + '.tsv'
        mn_ci_ouput_path = sample_dir + '/' + sample + mn_output_suffix + '_ci.tsv'
        repond_table_path = sample_dir + '/' + sample  + repond_table_suffix + '.tsv'
        repond_table_flat_path = sample_dir + '/' + sample + repond_table_suffix + '_flat.tsv'
        mn_graphml_ouput_path = sample_dir + '/' + sample + mn_output_suffix + '_mn.graphml'

        print('''
        Proceeding to spectral matching ...
        ''')

        spectrums_query = list(load_from_mgf(spectra_file_path))
        print('%s spectra were found in the query file.' % len(spectrums_query))

        with nostdout():
            spectrums_query = [metadata_processing(s) for s in spectrums_query]
            spectrums_query = [peak_processing(s) for s in spectrums_query]
        similarity_score = PrecursorMzMatch(tolerance=float(parent_mz_tol), tolerance_type="Dalton")
        chunks_query = [spectrums_query[x:x+1000] for x in range(0, len(spectrums_query), 1000)]

        for chunk in chunks_query:
            scores = calculate_scores(chunk, spectrums_db_cleaned, similarity_score)
            indices = np.where(np.asarray(scores.scores))
            idx_row, idx_col = indices
            scans_id_map = {}
            i = 0
            for s in chunk:
                scans_id_map[i] = int(s.metadata['scans'])
                i += 1
            if match_score == 'modifiedcosine':
                cosine = ModifiedCosine(tolerance=float(msms_mz_tol))
            else:
                cosine = CosineGreedy(tolerance=float(msms_mz_tol))
            data = []
            for (x,y) in tzip(idx_row,idx_col):
                if x<y:
                    msms_score, n_matches = cosine.pair(chunk[x], spectrums_db_cleaned[y])[()]
                    if (msms_score>float(min_cos)) & (n_matches>int(min_peaks)):
                        feature_id = scans_id_map[x]
                        data.append({'msms_score':msms_score,
                                    'matched_peaks':n_matches,
                                    'feature_id': feature_id,
                                    'reference_id':y + 1,
                                    'inchikey': spectrums_db_cleaned[y].get("name")})
            df = pd.DataFrame(data)
            df.to_csv(isdb_results_path, mode='a', header=not os.path.exists(isdb_results_path), sep = '\t')

        print('''
        Spectral matching done !
        ''')

        print('''
        Proceeding to Molecular Netorking computation...
        ''')      
        if mn_score == 'modifiedcosine':
            cosine = ModifiedCosine(tolerance=float(msms_mz_tol))
        else:
            cosine = CosineGreedy(tolerance=float(msms_mz_tol))

        scans_id_map[i] = int(s.metadata['scans'])
        #similarity_score = PrecursorMzMatch(tolerance=float(mn_parent_mz_tol), tolerance_type="Dalton")
        scores = calculate_scores(spectrums_query, spectrums_query, cosine)
        ms_network = SimilarityNetwork(identifier_key="scans", score_cutoff = mn_score_cutoff, top_n = mn_top_n, max_links = mn_max_links, link_method = 'mutual')
        ms_network.create_network(scores)
        ms_network.export_to_graphml(mn_graphml_ouput_path)
        components = connected_component_subgraphs(ms_network.graph)

        comp_dict = {idx: comp.nodes() for idx, comp in enumerate(components)}
        attr = {n: {'component_id' : comp_id} for comp_id, nodes in comp_dict.items() for n in nodes}
        comp = pd.DataFrame.from_dict(attr, orient = 'index')
        comp.reset_index(inplace = True)
        comp.rename(columns={'index': 'feature_id'}, inplace=True)
        count = comp.groupby('component_id').count()
        count['new_ci'] = np.where(count['feature_id'] > 1, count.index, -1)
        new_ci = pd.Series(count.new_ci.values,index=count.index).to_dict()
        comp['component_id'] = comp['component_id'].map(new_ci)
        spectrums_query_metadata_df = pd.DataFrame(s.metadata for s in spectrums_query)
        comp = comp.merge(spectrums_query_metadata_df[['feature_id', 'precursor_mz']], how='left')
        comp.to_csv(mn_ci_ouput_path, sep = '\t', index = False)

        print('''
        Molecular Networking done !
        ''')

        dt_isdb_results = pd.read_csv(isdb_results_path,
                                    sep='\t',
                                    usecols=['msms_score', 'feature_id', 'reference_id', 'inchikey'],
                                    error_bad_lines=False, low_memory=True)

        dt_isdb_results['libname'] = 'DNP_ISDB'

        dt_isdb_results.rename(columns={
            'inchikey': 'short_inchikey',
            'msms_score': 'score_input'}, inplace=True)

        clusterinfo_summary = pd.read_csv(mn_ci_ouput_path,
                                        sep='\t',
                                        usecols=['feature_id', 'precursor_mz', 'component_id'],
                                        error_bad_lines=False, low_memory=True)

        clusterinfo_summary.rename(columns={'precursor_mz': 'mz'}, inplace=True)

        cluster_count = clusterinfo_summary.drop_duplicates(
            subset=['feature_id', 'component_id']).groupby("component_id").count()
        cluster_count = cluster_count[['feature_id']].rename(
            columns={'feature_id': 'ci_count'}).reset_index()

        # ## we now merge this back with the isdb matched results 
        dt_isdb_results = pd.merge(dt_isdb_results, clusterinfo_summary, on='feature_id')

        db_metadata = pd.read_csv(metadata_path,
                                sep=',', error_bad_lines=False, low_memory=False)

        db_metadata['short_inchikey'] = db_metadata.structure_inchikey.str.split(
            "-", expand=True)[0]
        db_metadata.reset_index(inplace=True)

        print('Number of features: ' + str(len(clusterinfo_summary)))
        print('Number of MS2 annotation: ' + str(len(dt_isdb_results)))

        # Now we directly do the MS1 matching stage on the cluster_summary. No need to have MS2 annotations

        print('''
        Proceeding to MS1 annotation ...
        ''')
        super_df = []

        for i, row in tqdm(clusterinfo_summary.iterrows(), total=clusterinfo_summary.shape[0]):
            par_mass = clusterinfo_summary.loc[i, 'mz']
            df_0 = clusterinfo_summary.loc[[i], ['feature_id', 'mz', 'component_id']]
            df_1 = adducts_df[(adducts_df['min'] <= par_mass) & (adducts_df['max'] >= par_mass)]
            df_1['key'] = i
            df_1.drop(['min', 'max'], axis=1, inplace=True)
            df_tot = pd.merge(df_0, df_1, left_index=True, right_on='key', how='left')
            super_df.append(df_tot)

        df_MS1 = pd.concat(super_df, axis=0)
        del super_df

        df_MS1 = df_MS1.drop(['key'], axis=1).drop_duplicates(
            subset=['feature_id', 'adduct'])

        df_MS1['libname'] = 'MS1_match'

        print('''
        MS1 annotation done !
        ''')

        df_meta_2 = db_metadata[['short_inchikey', 'structure_exact_mass']]
        df_meta_2 = df_meta_2.dropna(subset=['structure_exact_mass'])
        df_meta_2 = df_meta_2.drop_duplicates(
            subset=['short_inchikey', 'structure_exact_mass'])

        df_meta_2 = df_meta_2.round({'structure_exact_mass': 5})
        df_MS1 = df_MS1.round({'exact_mass': 5})

        df_MS1_merge = pd.merge(df_MS1, df_meta_2, left_on='exact_mass',
                                right_on='structure_exact_mass', how='left')
        df_MS1_merge = df_MS1_merge.dropna(subset=['short_inchikey'])

        df_MS1_merge['match_mzerror_MS1'] = df_MS1_merge['mz'] - \
            df_MS1_merge['adduct_mass']
        df_MS1_merge = df_MS1_merge.round({'match_mzerror_MS1': 5}).astype({
            'match_mzerror_MS1': 'str'})

        df_MS1_merge = df_MS1_merge.drop(
            ['structure_exact_mass', 'adduct_mass', 'exact_mass'], axis=1)
        df_MS1_merge['score_input'] = 0

        # Merge MS1 results with MS2 annotations
        dt_isdb_results = pd.concat([dt_isdb_results, df_MS1_merge])

        print('Number of annotated features after MS1: ' +
            str(len(df_MS1_merge['feature_id'].unique())))

        print('Total number of MS1 and MSMS annotations: ' + str(len(dt_isdb_results)))

        # Rank annotations based on the spectral score

        dt_isdb_results["score_input"] = pd.to_numeric(
            dt_isdb_results["score_input"], downcast="float")
        dt_isdb_results['rank_spec'] = dt_isdb_results[['feature_id', 'score_input']].groupby(
            'feature_id')['score_input'].rank(method='dense', ascending=False)

        dt_isdb_results.reset_index(inplace=True, drop=True)

        # now we merge with the Occurences DB metadata after selection of our columns of interest

        cols_to_use = ['structure_inchikey', 'structure_inchi',
                    'structure_smiles', 'structure_molecular_formula',
                    'structure_exact_mass', 'short_inchikey', 'structure_taxonomy_npclassifier_01pathway', 
                    'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class',
                    'organism_name', 'organism_taxonomy_ottid',
                    'organism_taxonomy_01domain', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum',
                    'organism_taxonomy_04class', 'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_07tribe', 'organism_taxonomy_08genus', 'organism_taxonomy_09species', 'organism_taxonomy_10varietas' ]

        dt_isdb_results = pd.merge(
            left=dt_isdb_results, right=db_metadata[cols_to_use], left_on='short_inchikey', right_on='short_inchikey', how='outer')
        dt_isdb_results.dropna(subset=['feature_id'], inplace=True)

        print('Total number of annotations with unique Biosource/line: ' +
            str(len(dt_isdb_results)))


        # Now we want to get the taxonomic information for each of the samples
        # so we want to extract the species information from the metadata file
        samples_metadata[organism_header].dropna(inplace = True)
        samples_metadata[organism_header] = samples_metadata[organism_header].str.lower()
        species = samples_metadata[organism_header].unique()
        len_species = len(species)

        print("%s unique species have been selected from the metadata table." % len_species )

        species_tnrs_matched = OT.tnrs_match(species, context_name=None, do_approximate_matching=True, include_suppressed=False)

        with open(str(sample_dir +'species.json'), 'w') as out:
            sf = json.dumps(species_tnrs_matched.response_dict, indent=2, sort_keys=True)
            out.write('{}\n'.format(sf))

        with open(str(sample_dir +'species.json')) as tmpfile:
                jsondic = json.loads(tmpfile.read())

        json_normalize(jsondic)

        df_species_tnrs_matched = json_normalize(jsondic,
                    record_path=['results', 'matches']
                    )
        df_species_tnrs_unmatched = json_normalize(jsondic,
                    record_path=['unmatched_names']
                    )

        df_species_tnrs_matched.info()

        # We then want to match with the accepted name instead of the synonym in case both are present. 
        # We thus order by matched_name and then by is_synonym status prior to returning the first row.

        df_species_tnrs_matched.sort_values(['search_string', 'is_synonym'], axis = 0, inplace = True)
        df_species_tnrs_matched_unique = df_species_tnrs_matched.drop_duplicates('search_string', keep = 'first')

        # both df are finally merged
        merged_df = pd.merge(samples_metadata, df_species_tnrs_matched_unique, how='left', left_on=organism_header, right_on='search_string', indicator=True)

        # converting 'ott_ids' from float to int (check the astype('Int64') whic will work while the astype('int') won't see https://stackoverflow.com/a/54194908)
        merged_df['taxon.ott_id'] = merged_df['taxon.ott_id'].astype('Int64')

        # However, we then need to put them back to 
        merged_df['taxon.ott_id']
        ott_list = list(merged_df['taxon.ott_id'].dropna().astype('int'))

        taxon_info = []

        for i in ott_list:
            query = OT.taxon_info(i, include_lineage=True)
            taxon_info.append(query)

        tl = []

        for i in taxon_info:
            with open(str(sample_dir +'taxon_info.json'), 'w') as out:
                tl.append(i.response_dict)
                yo = json.dumps(tl)
                out.write('{}\n'.format(yo))

        with open(str(sample_dir +'taxon_info.json')) as tmpfile:
            jsondic = json.loads(tmpfile.read())

        df = json_normalize(jsondic)

        df_tax_lineage = json_normalize(jsondic,
                    record_path=['lineage'],
                    meta = ['ott_id', 'unique_name'],
                    record_prefix='sub_',
                    errors='ignore'
                    )

        # This keeps the last occurence of each ott_id / sub_rank grouping https://stackoverflow.com/a/41886945
        df_tax_lineage_filtered = df_tax_lineage.groupby(['ott_id', 'sub_rank'], as_index=False).last()

        #Here we pivot long to wide to get the taxonomy
        df_tax_lineage_filtered_flat = df_tax_lineage_filtered.pivot(index='ott_id', columns='sub_rank', values='sub_name')

        # Here we actually also want the lowertaxon (species usually) name
        df_tax_lineage_filtered_flat = pd.merge(df_tax_lineage_filtered_flat, df_tax_lineage_filtered[['ott_id', 'unique_name']], how='left', on='ott_id', )

        #Despite the left join ott_id are duplicated 
        df_tax_lineage_filtered_flat.drop_duplicates(subset = ['ott_id', 'unique_name'], inplace = True)

        # here we want to have these columns whatevere happens
        col_list = ['ott_id', 'domain', 'kingdom', 'phylum',
                                'class', 'order', 'family', 'tribe', 'genus', 'unique_name']

        df_tax_lineage_filtered_flat = df_tax_lineage_filtered_flat.reindex(columns=col_list, fill_value = np.NaN)

        # We now rename our columns of interest
        renaming_dict = {'domain': 'query_otol_domain',
                    'kingdom': 'query_otol_kingdom',
                    'phylum': 'query_otol_phylum',
                    'class': 'query_otol_class',
                    'order': 'query_otol_order',
                    'family': 'query_otol_family',
                    'tribe': 'query_otol_tribe',
                    'genus': 'query_otol_genus',
                    'unique_name': 'query_otol_species'}

        df_tax_lineage_filtered_flat.rename(columns=renaming_dict, inplace=True)

        # We select columns of interest 
        cols_to_keep = ['ott_id',
                    'query_otol_domain',
                    'query_otol_kingdom',
                    'query_otol_phylum',
                    'query_otol_class',
                    'query_otol_order',
                    'query_otol_family',
                    'query_otol_tribe',
                    'query_otol_genus',
                    'query_otol_species']

        df_tax_lineage_filtered_flat = df_tax_lineage_filtered_flat[cols_to_keep]

        # We merge this back with the samplemetadata only if we have an ott.id in the merged df 
        samples_metadata = pd.merge(merged_df[pd.notnull(merged_df['taxon.ott_id'])], df_tax_lineage_filtered_flat, how='left', left_on='taxon.ott_id', right_on='ott_id' )

        # Here we will add three columns (even for the simple repond this way it will be close to the multiple species repond)
        # these line will need to be defined as function arguments
        cols_att = ['query_otol_domain', 'query_otol_kingdom', 'query_otol_phylum', 'query_otol_class',
                    'query_otol_order', 'query_otol_family', 'query_otol_tribe', 'query_otol_genus', 'query_otol_species']
        for col in cols_att:
            dt_isdb_results[col] = samples_metadata[col][0]

        print('''
        Proceeding to taxonomically informed reponderation ...
        ''')

        cols_ref = ['organism_taxonomy_01domain', 'organism_taxonomy_02kingdom',  'organism_taxonomy_03phylum', 'organism_taxonomy_04class',
                    'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_07tribe', 'organism_taxonomy_08genus', 'organism_taxonomy_09species']

        cols_match = ['matched_domain', 'matched_kingdom', 'matched_phylum', 'matched_class',
                    'matched_order', 'matched_family', 'matched_tribe', 'matched_genus', 'matched_species']

        col_prev = None
        for col_ref, col_att, col_match in zip(cols_ref, cols_att, cols_match):
                dt_isdb_results[col_ref].fillna('Unknown', inplace=True)
                dt_isdb_results[col_att].fillna('Unknown', inplace=True)
                dt_isdb_results[col_ref] = dt_isdb_results[col_ref].apply(lambda x: [x])
                dt_isdb_results[col_att] = dt_isdb_results[col_att].apply(lambda x: [x])
                dt_isdb_results[col_match] = [list(set(a).intersection(set(b))) for a, b in zip(dt_isdb_results[col_ref], dt_isdb_results[col_att])] # Allows to compare 2 lists
                dt_isdb_results[col_match] = dt_isdb_results[col_match].apply(lambda y: np.nan if len(y)==0 else y)
                if col_prev != None:
                        dt_isdb_results[col_match].where(dt_isdb_results[col_prev].notnull(), np.nan)
                col_prev = col_match

        dt_isdb_results['score_taxo'] = dt_isdb_results[cols_match].count(axis=1)

        # Filter out MS1 annotations without a reweighting at the family level at least

        if polarity == 'pos':
            dt_isdb_results = dt_isdb_results[
                (dt_isdb_results['score_taxo'] >= min_score_ms1) | (
                dt_isdb_results['libname'] == 'ISDB')]
        else:
            dt_isdb_results = dt_isdb_results[
                (dt_isdb_results['score_taxo'] >= min_score_ms1) | (
                dt_isdb_results['libname'] == 'ISDB')]


        print('Total number of annotations after filtering MS1 annotations not reweighted at taxonomical level min: ' +
            str(len(dt_isdb_results)))

        print('Number of annotations reweighted at the domain level: ' +
            str(dt_isdb_results['matched_domain'].count()))
        print('Number of annotations reweighted at the kingom level: ' +
            str(dt_isdb_results['matched_kingdom'].count()))
        print('Number of annotations reweighted at the phylum level: ' +
            str(dt_isdb_results['matched_phylum'].count()))
        print('Number of annotations reweighted at the class level: ' +
            str(dt_isdb_results['matched_class'].count()))
        print('Number of annotations reweighted at the order level: ' +
            str(dt_isdb_results['matched_order'].count()))
        print('Number of annotations reweighted at the family level: ' +
            str(dt_isdb_results['matched_family'].count()))
        print('Number of annotations reweighted at the tribe level: ' +
            str(dt_isdb_results['matched_tribe'].count()))
        print('Number of annotations reweighted at the genus level: ' +
            str(dt_isdb_results['matched_genus'].count()))
        print('Number of annotations reweighted at the species level: ' +
            str(dt_isdb_results['matched_species'].count()))


        # we set the spectral score column as float
        dt_isdb_results["score_input"] = pd.to_numeric(
            dt_isdb_results["score_input"], downcast="float")
        # and we add it to the max txo score :
        dt_isdb_results['score_input_taxo'] = dt_isdb_results['score_taxo'] + \
            dt_isdb_results['score_input']


        dt_isdb_results['rank_spec_taxo'] = dt_isdb_results.groupby(
            'feature_id')['score_input_taxo'].rank(method='dense', ascending=False)

        dt_isdb_results = dt_isdb_results.groupby(["feature_id"]).apply(
            lambda x: x.sort_values(["rank_spec_taxo"], ascending=True)).reset_index(drop=True)

        # Get cluster Chemical class
        for col in ['structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class']:

            df = dt_isdb_results.copy()
            df = df.drop_duplicates(subset=['feature_id', col])
            df = df[df["component_id"] != -1]
            df = df[df.rank_spec_taxo <= top_N_chemical_consistency]
            df = df.groupby(
                ["component_id", col]
            ).agg({'feature_id': 'count',
                'rank_spec_taxo': 'mean'}
                ).reset_index(
            ).rename(columns={'feature_id': (col + '_count'),
                            'rank_spec_taxo': ('rank_' + col + '_mean')}
                    ).merge(cluster_count, on='component_id', how='left')

            df[('freq_' + col)] = df[(col + '_count')] / df['ci_count']
            df[(col + '_score')] = df[('freq_' + col)] / \
                (df[('rank_' + col + '_mean')]**(0.5))
            df = df.sort_values(
                (col + '_score'), ascending=False
            ).drop_duplicates(['component_id']
                            ).rename(columns={col: (col + '_consensus')})
            dt_isdb_results = dt_isdb_results.merge(
                df[[(col + '_consensus'), ('freq_' + col), 'component_id']], on='component_id', how='left')

        # Chemical consistency reweighting

        print('''
        Proceeding to chemically informed reponderation ...
        ''')

        dt_isdb_results['structure_taxonomy_npclassifier_01pathway_score'] = dt_isdb_results.apply(
            lambda x: 1 if x.structure_taxonomy_npclassifier_01pathway == x.structure_taxonomy_npclassifier_01pathway_consensus else 0, axis=1)
        dt_isdb_results['structure_taxonomy_npclassifier_02superclass_score'] = dt_isdb_results.apply(
            lambda x: 2 if x.structure_taxonomy_npclassifier_02superclass == x.structure_taxonomy_npclassifier_02superclass_consensus else 0, axis=1)
        dt_isdb_results['structure_taxonomy_npclassifier_03class_score'] = dt_isdb_results.apply(
            lambda x: 3 if x.structure_taxonomy_npclassifier_03class == x.structure_taxonomy_npclassifier_03class_consensus else 0, axis=1)

        dt_isdb_results['score_max_consistency'] = dt_isdb_results[[
            "structure_taxonomy_npclassifier_01pathway_score",
            "structure_taxonomy_npclassifier_02superclass_score",
            "structure_taxonomy_npclassifier_03class_score"
        ]].max(axis=1)

        dt_isdb_results['final_score'] = dt_isdb_results['score_input'] + dt_isdb_results['score_taxo'] + dt_isdb_results['score_max_consistency']

        dt_isdb_results['rank_final'] = dt_isdb_results.groupby(
            'feature_id')['final_score'].rank(method='dense', ascending=False)



        print('Number of annotations reweighted at the NPClassifier pathway level: ' +
            str(len(dt_isdb_results[(dt_isdb_results['structure_taxonomy_npclassifier_01pathway_score'] == 1)])))
        print('Number of annotations reweighted at the NPClassifier superclass level: ' +
            str(len(dt_isdb_results[(dt_isdb_results['structure_taxonomy_npclassifier_02superclass_score'] == 2)])))
        print('Number of annotations reweighted at the NPClassifier class level: ' +
            str(len(dt_isdb_results[(dt_isdb_results['structure_taxonomy_npclassifier_03class_score'] == 3)])))


        dt_isdb_results_chem_rew = dt_isdb_results.loc[(
            dt_isdb_results.rank_final <= int(top_to_output))]
        dt_isdb_results_chem_rew[["feature_id", "rank_final", "component_id"]] = dt_isdb_results_chem_rew[[
            "feature_id", "rank_final", "component_id"]].apply(pd.to_numeric, downcast='signed', axis=1)
        dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.sort_values(
            ["feature_id", "rank_final"], ascending=(False, True))
        dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.astype(str)


        # Here we would like to filter results when short IK are repeated for the same feature_id at the same final rank
        # see issue (https://gitlab.com/tima5/taxoscorer/-/issues/23)

        dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.drop_duplicates(subset=['feature_id', 'short_inchikey'], keep='first')

        dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.astype({'feature_id' : 'int64'})
       
        dt_isdb_results_chem_rew['lowest_matched_taxon'] = dt_isdb_results_chem_rew['matched_species']
        dt_isdb_results_chem_rew['lowest_matched_taxon'] = dt_isdb_results_chem_rew['lowest_matched_taxon'].replace('nan', np.NaN)
        col_matched = ['matched_genus', 'matched_tribe', 'matched_family', 'matched_order', 'matched_order', 'matched_phylum', 'matched_kingdom', 'matched_domain']
        for col in col_matched:
            dt_isdb_results_chem_rew[col] = dt_isdb_results_chem_rew[col].replace('nan', np.NaN)  
            dt_isdb_results_chem_rew['lowest_matched_taxon'].fillna(dt_isdb_results_chem_rew[col], inplace=True)

        annot_attr = ['rank_spec', 'score_input', 'libname', 'structure_inchikey', 'structure_inchi', 'structure_smiles', 'structure_molecular_formula', 'adduct',
                    'structure_exact_mass', 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class',
                    'query_otol_species', 'lowest_matched_taxon', 'score_taxo', 'score_max_consistency', 'final_score', 'rank_final']

        comp_attr = ['component_id', 'structure_taxonomy_npclassifier_01pathway_consensus', 'freq_structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass_consensus',
                    'freq_structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class_consensus', 'freq_structure_taxonomy_npclassifier_03class']


        col_to_keep = ['feature_id'] + comp_attr + annot_attr

        df4cyto_flat = dt_isdb_results_chem_rew[col_to_keep]

        all_columns = list(df4cyto_flat) # Creates list of all column headers
        df4cyto_flat[all_columns] = df4cyto_flat[all_columns].astype(str)

        gb_spec = {c: '|'.join for c in annot_attr}
        for c in comp_attr:
            gb_spec[c] = 'first'

        df4cyto = df4cyto_flat.groupby('feature_id').agg(gb_spec)

        df4cyto_flat.to_csv(repond_table_flat_path, sep='\t')

        df4cyto.to_csv(repond_table_path, sep='\t')

print('Finished in %s seconds.' % (time.time() - start_time))
