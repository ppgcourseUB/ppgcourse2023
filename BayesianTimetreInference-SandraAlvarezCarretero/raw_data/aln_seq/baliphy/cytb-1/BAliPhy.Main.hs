-- Use the program `brittany` (or `hindent`) to indent this file for readability
module Main where
import BAliPhy.ATModel
import BAliPhy.ATModel.DataPartition
import Bio.Alignment
import Bio.Alphabet
import IModel
import Parameters
import Probability
import SModel
import Tree
import qualified Data.Map as Map

sample_smodel_1 alpha = do {alpha_2 <- log_laplace 6.0 2.0
;kappa <- log_normal (log 2.0) 0.25
;pi <- symmetric_dirichlet_on (letters alpha) 1.0
;return ((hky85 alpha kappa (frequencies_from_dict alpha pi) & unit_mixture) & SModel.gamma_rates alpha_2 5, ["Rates.gamma:alpha" %=% alpha_2, "hky85:kappa" %=% kappa, "hky85:pi" %=% pi])
}

sample_imodel_1 = do {log_rate <- laplace (-4.0) 0.707
;mean_length <- shifted_exponential 10.0 1.0
;return (IModel.rs07 log_rate mean_length, ["rs07:log_rate" %=% log_rate, "rs07:mean_length" %=% mean_length])
}

sample_scale_1 = gamma 0.5 2.0

sample_branch_lengths tree = iid (numBranches tree) (gamma 0.5 (2.0 / intToDouble (numBranches tree)))

sample_topology_1 taxa = uniform_labelled_topology taxa

model taxa sequence_data imodel_training heat variable_alignment = do {topology1 <- sample_topology_1 taxa
;branch_lengths <- SamplingRate 0.0 (sample_branch_lengths topology1)
;scale1 <- sample_scale_1
;let {distances_part1 = listArray' (map (\x -> scale1 * x) branch_lengths)
}
;
;(smodel_1, log_smodel_1) <- sample_smodel_1 dna
;(_I1, log_I1) <- sample_imodel_1
;let {imodel_1 = _I1 topology1 heat imodel_training
}
;
;let {branch_hmms_part1 = branch_hmms imodel_1 distances_part1 13
}
;let {sequence_lengths_part1 = get_sequence_lengths (getAlphabet smodel_1) (sequence_data !! 0)
}
;alignment_on_tree_part1 <- random_alignment topology1 branch_hmms_part1 imodel_1 sequence_lengths_part1 variable_alignment
;
;return (ATModel topology1 [smodel_1] [scale1] branch_lengths Nothing [Partition smodel_1 (Just imodel_1) distances_part1 topology1 alignment_on_tree_part1 (Just branch_hmms_part1)], ["scale1" %=% scale1, "scale" %=% [scale1], "S1" %>% log_smodel_1, "I1" %>% log_I1])
}

main = do {let {imodel_training = Parameters.modifiable False
;heat = Parameters.modifiable 1.0
;variable_alignment = Parameters.modifiable True
;subst_root = Parameters.modifiable 13
}
;
;let {filenames = ["../unaln_raw.fasta"]
}
;let {filename_to_seqs = Map.fromList [(filename, load_sequences filename) | filename <- filenames]
}
;
;let {sequence_data1 = filename_to_seqs Map.! "../unaln_raw.fasta"
}
;
;let {sequence_data = [sequence_data1]
}
;
;let {taxa = map sequence_name sequence_data1
}
;
;(atmodel, loggers) <- random $ model taxa sequence_data imodel_training heat variable_alignment
;
;let {part1 = partitions atmodel !! 0
}
;let {(transition_ps_part1, cls_part1, ancestral_sequences_part1, likelihood_part1) = observe_partition_type_0 part1 sequence_data1 subst_root
}
;
;let {transition_ps = [transition_ps_part1]
}
;let {cond_likes = [cls_part1]
}
;let {anc_seqs = [ancestral_sequences_part1]
}
;let {likelihoods = [likelihood_part1]
}
;
;observe (fake_dist likelihoods) sequence_data
;
;return (ATModelExport atmodel transition_ps cond_likes anc_seqs likelihoods sequence_data imodel_training heat variable_alignment subst_root taxa, loggers)
}