import csv, os, time, sys, glob, numpy as np, scipy as sp, pandas as pd
from sklearn.model_selection import train_test_split, KFold
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.externals import joblib

# Set seed
np.random.seed(42)

def get_features(feat_fpath=None, thresh=1.00):
    df = pd.read_csv(feat_fpath, sep=',', header=None, index_col=0)
    data = df.transpose()
    return np.array(data), list(data.columns)


def pred_with_model(feat_fpath=None, out_pred_path=None, model_path='model.pkl'):
    itime = time.time()
    trained_regr = joblib.load(model_path)
    data, featnames = get_features(feat_fpath=feat_fpath)
    print(data.shape)
    print(str(time.time() - itime) + " sec to get data from file of computed features.")
    preds = trained_regr.predict(data)
    numcells_vec = data[:, featnames.index("NumCells")]
    if out_pred_path is not None:
        df = pd.DataFrame.from_dict({"NumCells": numcells_vec, "preds": preds})
        df.to_csv(out_pred_path, sep='\t', index=False, header=None)


if __name__ == "__main__":
    fargs = sys.argv[1:]
    print(fargs)
    pred_with_model(feat_fpath=fargs[0], out_pred_path=fargs[1], model_path='model.pkl')