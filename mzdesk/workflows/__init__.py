"""Multi-step workflows that run next to ESI MS/MS Annotation (mzdesk/annotation/pipeline.py).

common            the step runner and the generic result contract every workflow here returns
ei_ms_annotation  EI MS Annotation: one EI spectrum in, candidate metabolites out
targeted_method   Targeted Method: one compound in, a method sheet out
lipid_fragments   Lipid Fragment Explanation: a lipid name or structure + MS2 in, names, a structure and the explained peaks out

Each module exposes PIPELINE_ID, STEP_TITLES, Query, parse_inputs(inputs) -> Query and run_query(query, progress) -> result, which is what
mzdesk/api/app.py _PIPELINE_RUNNERS needs. Heavy engines are imported inside the steps, never at module level.
"""
