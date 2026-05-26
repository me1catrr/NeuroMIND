import numpy as np

from mne_brain.ica import compute_ica_features, evaluate_ica_components, suggest_rejected_components


def test_ica_feature_evaluation_shapes():
    rng = np.random.default_rng(42)
    component_maps = rng.normal(size=(4, 3))
    sources = rng.normal(size=(3, 1000))
    features = compute_ica_features(component_maps, sources, 500.0, ["Fz", "Cz", "T7", "Oz"])
    evaluated = evaluate_ica_components(features, artifact_thresh=10.0)

    assert features.shape[0] == 3
    assert "artifact_type" in evaluated.columns
    assert suggest_rejected_components(evaluated) == []

