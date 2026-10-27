"""Train-only spectral preprocessing inside a reproducible sklearn Pipeline.

Install the optional integration with: python -m pip install -e '.[ml]'
The tiny sequences here demonstrate the API; they are not a scientific dataset.
"""

import json

from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from protspect.sklearn import SequenceTransformer

train_sequences = ['AAAA', 'AAAC', 'AACA', 'CCCC', 'CCCA', 'CACC']
train_labels = [0, 0, 0, 1, 1, 1]
model = make_pipeline(
    SequenceTransformer(method='fourier', length=8, normalization='length',
                        standardize_channels=True, component=('amplitude', 'phase')),
    StandardScaler(),
    LogisticRegression(random_state=7),
)
model.fit(train_sequences, train_labels)
print('Predicted labels:', model.predict(['AAA', 'CCC']).tolist())
features = model.steps[0][1]
metadata = features.get_feature_metadata()
print('FFT length:', metadata['settings']['length'])
print('Feature count:', len(metadata['feature_names']))
# Store json.dumps(metadata) beside a feature table; persist the whole model
# with pickle/joblib to reuse its fitted grid, channel moments, and classifier.
assert json.loads(json.dumps(metadata)) == metadata
