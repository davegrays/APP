def fill_and_upsample(in_df, in_outcome, imp_meth, upsample):
    from fancyimpute import KNN
    from imblearn.over_sampling import SMOTE
    import pandas as pd
    import numpy as np

    if imp_meth == 'KNN':
        # Use 3 nearest rows which have a feature to fill in each row's missing features
        tmpdf = pd.DataFrame(np.round(KNN(k=3,verbose=False).complete(in_df)))
    elif imp_meth == 'median':
        imp = Imputer(strategy = "median")
        imp.fit(in_df)
        tmpdf = imp.transform(in_df)
    elif imp_meth == 'neg999':
        tmpdf = in_df.fillna(-999)
    else:
        tmpdf = in_df

    if upsample == True:
        # Use 5 nearest rows to upsample in order to have balanced classes
        sampler = SMOTE()
        X_res, y_res = sampler.fit_sample(tmpdf, in_outcome)
    else:
        X_res = tmpdf
        if isinstance(in_outcome, pd.Series):
            y_res = in_outcome.values
        else:
            y_res = in_outcome

    # convert predictors back to dataframe
    out_df = pd.DataFrame(X_res)
    out_outcome = np.round(y_res)
    if isinstance(in_df, pd.DataFrame):
        out_df.columns = in_df.columns
    return out_df, out_outcome