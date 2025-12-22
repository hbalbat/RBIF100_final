#!/usr/bin/env python3
# RBIF100 Week 7/8 - Final Project: Bioinformatics Data Analysis and Visualization
# Name: Helena Balbat
# making interactive version of a variant density histogram

# import packages
import pandas as pd
from dash import Dash, dcc, html
import plotly.express as px

# create dataframe with pandas of variants
df = pd.read_csv("data.csv")

# variant density histogram
fig_variant_density = px.histogram(
    df,
    x="position",
    nbins=50,  # number of bins
    title="APOE Variant Density Across Genomic Positions",
    labels={"position": "Genomic Position (bp)", "count": "Number of Variants"},
    color="variant_type",
    barmode="overlay"
)

# initialize dash app
app = Dash(__name__)

app.layout = html.Div([
    html.H1("Interactive APOE Variant Dashboard"),
    html.Div([
        html.H2("Variant Density Across Genomic Positions"),
        dcc.Graph(figure=fig_variant_density)
    ])
])

# run app
if __name__ == "__main__":
    app.run_server(debug=True)