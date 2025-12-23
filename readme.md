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
I thoroughly enjoyed working on this final project. It took some time, but I feel like it was a great way to apply everything that I've learned and practiced in this class. I would say the most rewarding part of the final project was creating my analysis section. I had the freedom to choose a gene of interest and analyze it in a way that made sense, which made the project engaging. I also felt rewarded when I debugged my code successfully. I feel like I understand RESTful APIs better, and I had less trouble generating a CSV file this time around. I wanted to take on the challenge of pushing all of my files to my pre-existing GitHub in a new repository, and that was where I had the biggest challenge. I have always worked with GitHub manually, aside from one or two times. I typically just uploaded files from my computer rather than pushing them from VS Code. This is where I consulted ChatGPT to help me figure that out a bit better. I ran into a few errors just in terms of syntax, so it helped me fix my commands in the terminal. I also am still new to Markdown files, so I had ChatGPT explain how to format my final summary report and this readme file (eg., bulleting format, how to make code blocks, how to italicize or bold things). I also used ChatGPT to help polish my final summary report, as I had written the entire report, but wanted to further polish the findings and interpretations sections. There were areas of blending, so I made sure to distinguish them better. For any time I used ChatGPT, there wasn't much need to "fact-check" it, because the information usually worked to do what I needed (markdown and pushing to GitHub). I made sure that it only provided me options to add/change in my final report as well, so I still had final say in what was written. In terms of the actual code, I felt comfortable enough referencing past assignments and notes to create my script.  

Working on this final project helped me see how computational tools tie together in bioinformatics research. Using scripting to retrieve data from NCBI via RESTful APIs showed me how researchers can access large datasets programmatically, rather than manually downloading them. I have almost always done the latter, so this is a strong skill to have and practice. Processing the data with Python and pandas gave me a better understanding of how to organize and clean biological datasets, while generating visualizations with matplotlib and seaborn made patterns in the data easier to interpret. I have worked with these libraries for multiple years now, but I always have things to learn and I am reinforcing the idea that these tools greatly benefit bioinformatics research/experimentation. Overall, this project proved how automation and coding skills enable reproducible and efficient analyses in today's world of bioinformatics.

If I had more time and resources to work on this analysis, I would most likely generate more figures or even explore more ways to manipulate/analyze my dataset. I unfortunately ran into a common issue where NCBI variants report being SNPs even if the reference allele and alternate allele appear to be the same in the requested data. I found that this happens in dbSNP because the SNP database aggregates reported variants across populations and genome assemblies. A variant may be listed in some individuals even if it matches the reference genome in others. I found that this is typical when working with public variant databases. Although I would have to re-analyze results, I am curious as to how my data would change if I selected a new gene of interest. I made sure to create reusable functions for the analysis, so I would only have to adjust titles for new gene of interest (and perhaps a few other things in the visualizations section). If I had been able to create a larger dataset, there may have been more SNPs without the same reference and alternate alleles, or some with insertions as well. I already created a larger dataset than required, but I am curious as to what would happen if I expanded the dataset and refined it to those without the same alleles. I'd also be interested in practicing with Plotly and Dash more to create more interactive visualizations.