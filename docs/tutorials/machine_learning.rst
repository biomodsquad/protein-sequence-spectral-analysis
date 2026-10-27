Train within a reproducible Pipeline
====================================

Put ``SequenceTransformer`` inside the sklearn Pipeline. It learns the optional
channel means/scales and implicit grid from the training fold only. The feature
scaler and classifier must also fit within that fold. A fixed explicit length
is useful when you can specify sequence coverage independently of test outcomes.

A small runnable example
------------------------

These hand-constructed classes demonstrate the API, not antimicrobial activity.

.. testcode::

   import json
   import numpy as np
   from sklearn.model_selection import GridSearchCV, StratifiedKFold
   from sklearn.pipeline import Pipeline
   from sklearn.preprocessing import StandardScaler
   from sklearn.linear_model import LogisticRegression
   from protspect.sklearn import SequenceTransformer
   from protspect.embeddings import kidera_factors

   train = ["KKAA", "KAKA", "AKKK", "RRKA", "KRAA", "ARKK",
            "DDAA", "DADA", "ADDD", "EEDA", "DEAA", "AEDD"]
   y = np.array([1] * 6 + [0] * 6)
   pipeline = Pipeline([
       ("features", SequenceTransformer(embedding=kidera_factors(), method="fourier",
           length=16, normalization="length", standardize_channels=True)),
       ("scale", StandardScaler()),
       ("model", LogisticRegression(max_iter=1000, random_state=7)),
   ])
   folds = StratifiedKFold(3, shuffle=True, random_state=7)
   search = GridSearchCV(pipeline, {"model__C": [.1, 1.]}, cv=folds,
                          scoring="balanced_accuracy", error_score="raise")
   search.fit(train, y)
   fitted = search.best_estimator_.named_steps["features"]
   before = fitted.get_feature_metadata()
   probabilities = search.predict_proba(["KKRA", "DDEA"])[:, 1]
   assert probabilities.shape == (2,) and np.isfinite(probabilities).all()
   assert fitted.get_feature_metadata() == before
   assert json.loads(json.dumps(before)) == before

Tuning transforms
-----------------

Use a list of parameter dictionaries to associate methods with their own options.
Replace the grid above with a controlled, equal-budget comparison:

.. code-block:: python

   grid = [
       {"features__method": ["fourier"], "model__C": [.1, 1.]},
       {"features__method": ["z"],
        "features__transform_options": [{"radius": 1.05}], "model__C": [.1, 1.]},
       {"features__method": ["mellin"],
        "features__transform_options": [{"mellin_c": 0.0}], "model__C": [.1, 1.]},
       {"features__method": ["wavelet"],
        "features__transform_options": [{"wavelet": "haar", "level": 2,
                                          "mode": "periodization"}], "model__C": [.1, 1.]},
   ]

Notebook 03 in :doc:`../examples/index` runs a four-transform search on synthetic
sequences and scores only the selected model on its held-out synthetic set.
For a real MLP comparison use :doc:`antimicrobial_mlp`.

Freeze and reuse
----------------

``get_feature_metadata`` returns independent JSON-safe embedding scores, grids,
preprocessing, moments, selection, and feature names. Persist the whole fitted
Pipeline with pickle/joblib for prediction, not just this metadata. Load serialized
models only from trusted sources and record the library environment. Re-fitting
is required to apply constructor changes; ``set_params`` invalidates fitted state.

A held-out chain longer than the fitted grid raises instead of changing columns.
Do not fix this by inspecting test scores or silently truncating. Specify a larger
grid from the application requirements and refit the whole training workflow.
Use family-aware splits or external data to assess biological generalization.
