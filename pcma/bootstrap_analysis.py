from joblib import Parallel, delayed
import statsmodels.api as sm


def bootstrap_p_values(media_df, encoded_y, x_col, mediator_col, n_bootstrap):

    def bootstrap_iteration(media_df, encoded_y, x_col, mediator_col):
        combined_df = media_df.copy()
        combined_df['encoded_y'] = encoded_y
        sampled_df = combined_df.sample(frac=1, replace=True)

        sampled_encoded_y = sampled_df['encoded_y']
        sampled_df = sampled_df.drop('encoded_y', axis=1)

        model2_bootstrap = sm.OLS(sampled_df[mediator_col],
                                  sm.add_constant(sampled_df[x_col])).fit()
        coef_a_bootstrap = model2_bootstrap.params.iloc[1]

        model3_bootstrap = sm.Logit(
            sampled_encoded_y,
            sm.add_constant(sampled_df[[x_col, mediator_col]])).fit(disp=False)
        coef_b_bootstrap = model3_bootstrap.params.iloc[2]

        t_stat_a_bootstrap = coef_a_bootstrap / model2_bootstrap.bse.iloc[1]
        t_stat_b_bootstrap = coef_b_bootstrap / model3_bootstrap.bse.iloc[2]

        return t_stat_a_bootstrap * t_stat_b_bootstrap

    bootstrap_z_values = Parallel(n_jobs=-1)(
        delayed(bootstrap_iteration)(media_df, encoded_y, x_col, mediator_col)
        for _ in range(n_bootstrap))

    return bootstrap_z_values
