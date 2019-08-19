# Sampler
# In: parameter space
# Out: parameter sets
#
#
# For step in budget:
#
#   For rep in reps:
#
#       Model
#       In: parameter set[step]
#       Out: model data
#
#       Evaluation
#       In: model data & either true data or stylised facts
#       Out: labelled data (for one run)
#
#   Result averaging
#   In: Labels for each run
#   Out: Aggregated labels (one for each criterion)
#
#   If step is big enough to run surrogate with:
#
#       Surrogate
#       In: labelled data
#       Out: predictions for unlabelled data
#
#       Selector
#       In: predictions for unlabelled data
#       Out: unlabelled parameter set
