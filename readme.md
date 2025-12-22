# RBIF100 Weeks 7/8 - Final Project: Bioinformatics Data Analysis and Visualization
## Name: Helena Balbat
### Note: All files have been pushed to GitHub

See GitHub link: https://github.com/hbalbat/RBIF100_final 

## Purpose:
This final assignment is intended to apply all the skills we've learned throughout this course to analyze a real-world biological dataset. This includes scripting, automation, APIs, bioinformatics data processing, and visualization. These Python scripts use libraries like matplotlib, pandas, and NumPy to analyze variant types of the *Homo sapiens* apolipoprotein E (APOE) gene. The final_analysis.py script generates a data CSV file and three figures. I also took the challenge of pushing all files to a GitHub repository (using my pre-existing account), as well as practicing using Plotly and Dash. 


## Included Files:
* final_analysis.py : Python3 script that includes steps 1-3 of the final assignment. Generates data.csv
* interactive_visualizations.py : Python3 script that utilizes Plotly and Dash to generate an interactive histogram of variant density data across the apolipoprotein E (APOE) gene
* analysis_report.md : final summary report of data analysis. Includes three generated figures 

## Generated Files
* data.csv : generated CSV file of data collected with a RESTful API request from NCBI. Includes information on variants of the *Homo sapiens* APOE gene
* figure1.png : bar chart of variant type distribution across the APOE gene
* figure2.png : bar chart of variant density distribution across the APOE gene
* figure3.png : boxplot of variant density distribution across APOE genomic bins


## How to Run:
In remote class server, type the following:
        
    cd /home/balbh/week7
    python3 final_analysis.py

This will output general print statements and results of the performed analyses. The script will output necessary data to data.csv

To generate the interactive histogram, type the following:
        
    pip install dash plotly
    python3 interactive_visualizations.py

When the pop-up appears on your screen, you can open the visualization in a browser.


## Reflection and AI Use Disclosure:
Which part of this final project did you find most interesting or rewarding? Which was the most challenging?
How did integrating scripting, API use, data handling, and visualization impact your understanding of bioinformatics research?
How would you extend or improve your analysis further with more time or resources?