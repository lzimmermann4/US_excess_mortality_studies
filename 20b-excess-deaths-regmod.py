import numpy as np
import pandas as pd
from regmod.model import PoissonModel
from regmod.data import Data
from regmod.parameter import Parameter
from regmod.variable import Variable
from regmod.function import SmoothFunction
from patsy import dmatrix
import matplotlib.pyplot as plt
import os
import random
import statsmodels.api as sm
from sklearn.metrics import mean_squared_error

### This script was adapted from 02b_fit_regmod_model.py in IHME's GitHub repository, in order to perform the four base models which rely on regmod.
data_wk = pd.read_csv("[insert file path to working directory]/data/data_all_cause/USA.csv")

config = {
    'n_draws': 1000,  
    'summary_probs': [0.025, 0.975]  
}

# Convert weekly data to monthly data
data_wk['time_start'] = pd.to_datetime(data_wk['time_start'])
data_wk['month'] = data_wk['time_start'].dt.month

data = data_wk.groupby(['year', 'month']).agg({
    'deaths': 'sum',
    'population': 'max'
}).reset_index()

data['time'] = np.arange(1,len(data)+1)
data['mortality_rate'] = data['deaths'] / data['population']
data['offset'] = np.log(data['population'])

# Tail sizes in months (6, 12, 18, 24 months)
tail_sizes = [6, 12, 18, 24]

# Directory to save results
output_dir = "regmod_results"
os.makedirs(output_dir, exist_ok=True)

# Initial split into reference and prediction periods
ref_data = data[data['year'] <= 2019].copy()
pred_data = data[data['year'] >= 2020].copy()

# Function to generate knots with a given number per year
# note: use 15 knots across period to strike a balance 
#       between 3 knots * 6 years = 18 knots and 3 * 4 years= 12 knots 
#       to have ref and pred periods have the same number of knots
def get_time_knots_adj(df, tail_size, num_knots=15):
    # Place knots evenly across the reference period with consideration of tail
    time_min = df['time'].min()
    time_max = df['time'].max() - tail_size
    
    boundary_knots = [time_min, time_max + tail_size]
    internal_knots = np.linspace(time_min, time_max, num_knots - 2)
    knots = np.concatenate([boundary_knots[:1], internal_knots, boundary_knots[1:]])
    return knots

# note: the IHME's function for knots does not work when using reference and prediction periods of difference sizes
#       however, they appear to be using about 3 knots per year (2.75) and so we replicate this above
def get_time_knots(df, tail_size, units_per_year=12, knots_per_year=0.5):
    time_min = df['time'].min()
    time_max = df['time'].max()

    body_size = time_max - time_min - tail_size + 1 
    num_body_knots = int(knots_per_year*body_size/units_per_year) + 1
    if num_body_knots < 2:
        time_knots = np.array([time_min, time_max])
    else:
        time_knots = np.hstack([
            np.linspace(time_min, time_max - tail_size, num_body_knots),
            time_max
        ])
    return time_knots

for tail_size in tail_sizes:
    print(f"Running model with a tail size of {tail_size}")

    knots = get_time_knots_adj(df=ref_data, tail_size=tail_size)
    spline_matrix = dmatrix("bs(time, knots=knots[1:-1], degree=3, include_intercept=False)",
                            {"time": ref_data['time']},
                            return_type='dataframe')

    for col in spline_matrix.columns:
        ref_data[col] = spline_matrix[col]
    print(knots)

    # Prepare data using regmod's Data class
    regmod_data = Data(
        df=ref_data,
        col_obs='mortality_rate',
        col_covs=spline_matrix.columns.tolist(),
        col_offset='offset'
    )

    variables = [Variable(name=col) for col in spline_matrix.columns]

    # Define smooth exponential function 
    exp_function = SmoothFunction(
        name="exp", 
        fun=lambda x: np.exp(x),  
        inv_fun=lambda x: np.log(x),  
        dfun=lambda x: np.exp(x),      
        d2fun=lambda x: np.exp(x)      
    )

    # Define parameters with regmod Parameter class
    params = [
        Parameter(
            name="lam", 
            variables=variables,
            inv_link=exp_function 
        )
    ]

    # Initialize Poisson model with parameters
    model = PoissonModel(data=regmod_data, params=params)
    model.fit()

    knots = get_time_knots_adj(df=pred_data, tail_size=tail_size)
    print(knots)
    spline_matrix_pred = dmatrix("bs(time, knots=knots[1:-1], degree=3, include_intercept=False)",
                                       {"time": pred_data['time']},
                                       return_type='dataframe')
    
    for col in spline_matrix_pred.columns:
        pred_data[col] = spline_matrix_pred[col]

    predicted_lambda = np.exp(np.dot(spline_matrix_pred.values, model.opt_coefs) + pred_data['offset'])

    # Draw samples from Poisson distribution with lambda
    predicted_draws = np.random.poisson(lam=predicted_lambda[:, np.newaxis], size=(len(pred_data), config['n_draws']))

    # Calculate excess deaths for each draw
    excess_deaths_draws = pred_data['deaths'].values[:, np.newaxis] - predicted_draws

    pred_data_long = pd.DataFrame({
        'year': np.repeat(pred_data['year'].values, config['n_draws']),
        'month': np.repeat(pred_data['month'].values, config['n_draws']),
        'time': np.repeat(pred_data['time'].values, config['n_draws']),
        'deaths': np.repeat(pred_data['deaths'].values, config['n_draws']),
        'population': np.repeat(pred_data['population'].values, config['n_draws']),
        'draw': np.tile(np.arange(1, config['n_draws'] + 1), len(pred_data)),
        'expected_deaths': predicted_draws.flatten(),
        'excess_deaths': excess_deaths_draws.flatten(),
        'tail_size': tail_size
    })
    pred_data_long['model_type'] = f"regmod {tail_size}"
    pred_data_long.head(5)

    # Save the results to a CSV file
    output_file = os.path.join(output_dir, f"regmod_ts_{tail_size}_months.csv")
    pred_data_long.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

    ## Prepare data for plots
    plot_dat = pred_data_long.groupby(['year', 'month', 'time']).agg({
        'deaths': 'mean',
        'expected_deaths': 'mean',
        'excess_deaths': 'mean'
    }).reset_index()
    plot_dat['year_month'] = plot_dat['year'].astype(str) + '-' + plot_dat['month'].astype(str).str.zfill(2)
    
    ## Examine plots
    # Plot actual vs. expected deaths
    plt.figure(figsize=(10, 6))
    plt.plot(plot_dat['year_month'], plot_dat['deaths'], label='Observed Deaths', marker='o', color='darkgrey')
    plt.plot(plot_dat['year_month'], plot_dat['expected_deaths'], label='Expected Deaths', linestyle='--', color='orange')
    plt.xlabel('Year-Month')
    plt.ylabel('Number of Deaths')
    plt.xticks(rotation=65)
    plt.legend()
    plt.title(f'Piece-wise Poisson Regression with Tail Size of {tail_size} Months')
    plt.tight_layout()

    # Save plot as PNG
    plot_file = os.path.join(output_dir, f"expected_deaths_ts_{tail_size}_months.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()
    print(f"Plot saved to {plot_file}")

    # Plot excess deaths
    plt.figure(figsize=(10, 6))
    plt.plot(plot_dat['year_month'], plot_dat['excess_deaths'], label='Excess Deaths', marker='o')
    loess = sm.nonparametric.lowess
    smoothed_excess_deaths = loess(plot_dat['excess_deaths'], plot_dat['time'], frac=0.3)
    plt.plot(plot_dat['year_month'], smoothed_excess_deaths[:, 1], label='LOESS Smoothed Excess Deaths', linestyle='-', color='lightblue')
    plt.xlabel('Year-Month')
    plt.ylabel('Excess Deaths')
    plt.xticks(rotation=65)
    plt.legend()
    plt.title(f'Excess Deaths with Tail Size of {tail_size} Months')

    # Save plot as PNG
    plt.tight_layout()
    plot_file = os.path.join(output_dir, f"excess_deaths_ts_{tail_size}_months.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()
    print(f"Plot saved to {plot_file}")
