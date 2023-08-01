NUPACK_WRAP_ENUM(Ensemble);
NUPACK_WRAP_ENUM(Parameters);

NUPACK_WRAP(HeartbeatProgress, complete, frontend);
NUPACK_WRAP(BatchProgress, current, total);
NUPACK_WRAP(DefectProgress, defect, weighted_defect, new_upload, root, trial);
NUPACK_WRAP(ExceptionProgress, exception);
NUPACK_WRAP_VARIANT(ProgressContent, heartbeat, batch, defect, exception);
NUPACK_WRAP(Progress, id, time, content);

NUPACK_WRAP(Strand, index, name, sequence);
NUPACK_WRAP(Model, kelvin, sodium, magnesium, parameters, ensemble, filename);
NUPACK_WRAP(Complex, strand_index);

NUPACK_WRAP(PairProbability, diagonal, rows, cols, values);
NUPACK_WRAP(ComplexSet, strands, complex_indices);
NUPACK_WRAP(PairsJob, model, complexes);
NUPACK_WRAP(MFEJob, model, complexes);

NUPACK_WRAP(PairsResult, pairs, log_partition_function);
NUPACK_WRAP(MFEResult, mfe_structure, mfe);

NUPACK_WRAP(PairsResults, results);
NUPACK_WRAP(MFEResults, results);

NUPACK_WRAP(DesignDomain, name, sequence);
NUPACK_WRAP(DesignStrand, name, domains);
NUPACK_WRAP(DesignComplex, name, strands, structure);
NUPACK_WRAP(DesignTube, name, complexes, concentrations, includes, excludes, max_size);
NUPACK_WRAP(DesignOptions, objective_weight, seed, f_stop, f_passive, h_split, n_split, f_split,
     f_stringent, dg_clamp, m_bad, m_reseed, m_reopt, f_redecomp, f_refocus, f_sparse, checkpoint_time, max_time, defect_time, wobble_mutations);

NUPACK_WRAP(MatchConstraint, left, right);
NUPACK_WRAP(ComplementarityConstraint, left, right, wobble_mutations);
NUPACK_WRAP(DiversityConstraint, domains, word_length, min_nucleotide_types);
NUPACK_WRAP(WindowConstraint, domains, sequences);
NUPACK_WRAP(SequenceList, sequences);
NUPACK_WRAP(LibraryConstraint, domains, libraries);
NUPACK_WRAP(PatternConstraint, domains, patterns);
NUPACK_WRAP(SimilarityConstraint, domains, reference, min, max);
NUPACK_WRAP(SsmConstraint, complexes, word_size);
NUPACK_WRAP(EnergyMatchConstraint, domains, energy_ref);
NUPACK_WRAP_VARIANT(HardConstraint, match, complementarity, pattern, diversity, library, window, similarity);
NUPACK_WRAP_VARIANT(ConstraintObjective, pattern, similarity, ssm, energy_match);
NUPACK_WRAP(SoftConstraint, constraint, weight);

NUPACK_WRAP(DesignWeight, tube, complex, strand, domain, weight);
NUPACK_WRAP(DesignJob, model, domains, strands, complexes, tubes, options, hard_constraints, soft_constraints, weights, trial);

NUPACK_WRAP(DesignComplexResult, name, log_partition_function, defect, normalized_defect, pair_probabilities, energy);
NUPACK_WRAP(DesignTubeComplexResult, name, concentration, target_concentration, defect, structural_defect, concentration_defect, normalized_defect_contribution);
NUPACK_WRAP(DesignTubeResult, name, nucleotide_concentration, defect, normalized_defect, complexes);

NUPACK_WRAP(DesignPartition, mask, deflate);
NUPACK_WRAP(DesignStats, num_leaf_evaluations, num_reseeds, num_redecompositions, offtargets_added_per_refocus, design_time, analysis_time, final_psi, seed);
NUPACK_WRAP(DefectHistory, defects, wall_times, root);
NUPACK_WRAP(DesignResult, domains, unweighted_defects, weighted_defects, complexes, tubes, stats, history, trial);

NUPACK_WRAP(EnergyJob, strands, model, complex, structure);

NUPACK_WRAP(ExceptionResult, exception);
NUPACK_WRAP(BackendCancelledResult, timestamp);

NUPACK_WRAP_VARIANT(BackendJobInput, pairs_job, mfe_job, design_job);
NUPACK_WRAP(BackendJob, id, frontend_job_id, slots, input);

NUPACK_WRAP_VARIANT(BackendJobOutput, pairs_results, mfe_results, design_results, exception, cancelled);
NUPACK_WRAP(BackendResult, id, resources_used, output);

NUPACK_WRAP_VARIANT(CheckpointContent, design_checkpoint);
NUPACK_WRAP(Checkpoint, time, id, content);

/*************************************************************************/

NUPACK_WRAP(FrontendAnalysisTube, log_concentrations, strands, name, complexes, includes, excludes);
NUPACK_WRAP(FrontendAnalysisJob, tubes, models);
NUPACK_WRAP(FrontendConcentrationJob, tubes, models, base_id);
NUPACK_WRAP(FrontendUtilityPairsJob, model, strands, complex, structure);
NUPACK_WRAP(FrontendUtilityDesignJob, model, domains, strands, complex, options);
NUPACK_WRAP_VARIANT(FrontendUtilityJob, design, analysis); // id, parent_id
NUPACK_WRAP(FrontendDesignJob, spec, trials);
NUPACK_WRAP_VARIANT(FrontendJob, analysis, concentration, design, utility); // id, parent_id
NUPACK_WRAP(FrontendAnalysisConditions, tube, model);
NUPACK_WRAP(FrontendAnalysisTubeResult, conditions, complexes, pairs, mfe);
NUPACK_WRAP(FrontendAnalysisResults, tube_result);
NUPACK_WRAP(FrontendCancellationRequest, id, time, audit_reason);
NUPACK_WRAP(FrontendCancellationResult, request); // missing cancel time on purpose
NUPACK_WRAP(FrontendDesignResults, trial_names, results); // stop_request
NUPACK_WRAP(FrontendUtilityResult, structure, energy, pairs);
NUPACK_WRAP_VARIANT(FrontendResult, analysis, design, exception, cancellation, utility);
