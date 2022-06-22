import ctypes
import xgboost
import math
import pandas as pd
import numpy as np
import json
import pickle
import types

def comp_model(A, B):
    def _comp_model_var(k, v, B):
        if type(v) == ctypes.c_void_p:
            return True # Assume all C pointers are True
        elif type(v) == np.ndarray:
            return np.array_equal(v, vars(B)[k])
        else:
            return v == vars(B)[k]


    def xgb_to_dict(model):
         with tempfile.NamedTemporaryFile(suffix=".json") as fp:
            model.save_model(fp.name)
            with open(fp.name, "r") as f:
                return json.load(f)

    if type(A) == xgboost.sklearn.XGBClassifier:
        if type(B) != xgboost.sklearn.XGBClassifier:
            return False
        return {"same_xgb": xgb_to_dict(A) == xgb_to_dict(B)}

    return {k: _comp_model_var(k, v, B) for k, v in vars(A).items()}

def test_models(A,B):
    def _comp_models(A, B):
        same_type = all([type(x) == type(y) for x, y in zip(A, B)])
        if not same_type:
            return False
        same_models = pd.DataFrame([comp_model(x, y) for x, y in zip(A, B)])
        return all(same_models)
    
    def _comp_var(k, v, B):
        if type(v) in [int, float, str, bool, types.FunctionType]:
            return v == vars(B)[k]
        elif v is None or vars(B)[k] is None:
            return v is None and vars(B)[k] is None
        elif k == "models":
            return _comp_models(v, vars(B)[k])
        elif k == "model":
            same_type = type(v) == type(vars(B)[k])
            if not same_type:
                return False
            elif type(v) == xgboost.sklearn.XGBClassifier:
                return all(comp_model(v, vars(B)[k]))
            # CRF and CNN not implemented and will return None
    
    return {k: _comp_var(k, v, B) for k, v in vars(A).items()
            if k not in ['n_jobs', 'time']}

pkl_sema4 = 'BioMe_array-sema4_GDA_chr17/models/model_chm_17/model_chm_17.pkl'
pkl_nyscf = 'NYSCF_array_chr17/models/model_chm_17/model_chm_17.pkl'

with open('output/' + pkl_sema4, 'rb') as f:
    Sema4 = pickle.load(f)
with open('output/' + pkl_nyscf, 'rb') as f:
    NYSCF = pickle.load(f)

def test_all(A, B):
    def cm_eq(x):
        return np.array_equal(x[0][0], x[1][0]) and x[1][1] == x[1][1]
    def _comp_var(k, v, B):
        if k == "Confusion_Matrices":
            AB = pd.DataFrame([v, vars(B)[k]]).to_dict(orient="list")
            return {kk: cm_eq(vv) for kk, vv in AB.items()}
        elif type(v) in [int, float, str, bool, dict, types.FunctionType]:
            return v == vars(B)[k]
        elif v is None or vars(B)[k] is None:
            return v is None and vars(B)[k] is None
        elif k in ["base", "smooth"]:
            return test_models(v, vars(B)[k])
        elif type(v) == np.ndarray:
            return np.array_equal(v, vars(B)[k])
        elif type(v) == pd.core.frame.DataFrame:
            return v.equals(vars(B)[k])
    
    return {k: _comp_var(k, v, B) for k, v in vars(A).items()
            if k not in ['n_jobs', 'time']}

def rm_true_nested(x):
    y = x.copy()
    for k, v in y.copy().items():
        if v == True:
            del y[k]
        elif type(v) == dict:
            nret = test_nested(v)
            if len(nret) == 0:
                del y[k]
            else:
                y[k] = nret
    return y


all_results = test_all(NYSCF, Sema4)

issues = rm_true_nested(all_results)

if len(issues) == 0:
    print("Models are identical")
else:
    print("The following model components differ")
    print(issues)