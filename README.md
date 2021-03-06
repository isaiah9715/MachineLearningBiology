# Check out the "Updated_Analysis" file. Here, I use LASSO regularization to predict cancer stages and reduce the dimisionality of the model from 18,427 genes to 133 genes! 

# Stage Dependent Modeling of Cervical Cancer from Gene Expression Data (Machine Learning on Biological Data using R)

## How this README is organized (If you're only interested in the code it is TheCode.md file)
1. Project Description
2. Project Inspiration 
3. A little bit about LASSO
4. The Workflow of the project
5. Notes
6. Description of Code and Data
8. Future suggestions/goals 

## Project Description:
I'm being funded to help the [Wallace Lab](https://wallacelabksu.weebly.com/) at Kansas State University. Specifically, I'm applying [LASSO](https://www.jstor.org/stable/2346178?seq=1), which is an efficient machine learning model used for feature selection and regularization. I also did boxplot analysis of the genes, and then I made a heatmap of the data by using complete linkage. 

__The main goal of this project is to create a predictive model that will tell us which stage of cancer a tissue sample is in based off RNA sequenced data from genes__

## Project Inspiration
Prophylactic vaccines can prevent the human papillomavirus (HPV) infections that cause cervical cancers (CaCx). Yet, CaCx still kills a woman every 90 seconds. Vaccine hesitancy hinders the herd immunity and broad protection in wealthy countries. Some of the hesitancy comes from the very low rates of CaCx thanks to widespread frequent pap smear screening that are quite effective at providing early interventions that prevent these malignancies. Indeed, the 14 million new HPV infections each year in the US, only result in 13,170 cases of CaCx. In developing countries, where the majority of CaCx deaths occur, fiscal and infrastructural constraints further any hope for a successful vaccination campaign and present significant barriers to pap smear screening. 
While these screening modalities are quite effective at detecting pre-malignant lesions, they could be done less frequently if we could identify women who are at higher risk of progression. Most HPV infections and many premalignant lesions would be cleared by the host immune response without ever posing a cancer risk. However, since there are no established biomarkers for triaging women based on the relative likelihood of progression, paps smears have become routine.  If women at greater risk could be identified, paps smear screening could be optimized by less frequent screening of low risk cohorts. In higher socioeconomic groups, this would decrease the stress and physical discomfort associated with the exam. The potential rewards in low and middle-income nations would be far greater. In these communities, reducing the number of interventions would allow limited resources to be better allocated to people with the greatest need.

### LASSO
I'm using Multinomial LASSO which is found in the GLMNET package description [here](https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html). It is specifically under the Multinomial Models section of the webpage. 

### Why LASSO?
LASSO is great for parsimonious modeling when dealing with high-dimesional n<<p datasets or just multi-dimesional datasets because of the shrinkage it conducts (feature selection). This helps us deal with overfitting. In my case, the LASSO model is doing feature selection only since n>p. 

## Workflow
I've been working on this project for almost 6 months now. I have to check in with my research team every week to show them what I'm working on since I am being funded. This also helps me a great deal because they let me know whether my analysis makes sense, because I don't know too much about the biology happeneing. __I'm just the statistician in this research team, which is an honor in itself.__

This project started off with gathering the data needed to conduct any analysis. We screened and filtered through hundreds of thousands of observations coming from __7 different NCBI GEO Cervical Cancer Databases.__ 

After that, the team had an ample supply of data to conduct analysis but I also had to __clean the data__. First, I normalized and merged all the data. Since not all of the observation measurments were derived the same way normalizing was necessary. Then I went through the entire dataset and let the team know if I seen anything weird/concerning: Lot's of NA's, 0's, Duplicates, etc. 

The dataset we obtained had a little more than 25,000 columns and 400 rows. This is a high-dimesional data set. The team then told me they were concerned with a specific set of genes, 23 genes to be  exact. I don't know why they chose the ones they did, but they know more about the biology than me. So, after all the data preparation the dataset being used for this code had 22 columns and over 260 rows. 22 columns was supposed to be 23, but the data for some of the genes wasn't suitable for analysis or just wasn't there. 

Finally, we had the dataset we wanted to work with. I made some heatmaps, boxplots, and applied LASSO to the dataset. Most of the time was spent getting my data ready vs actually doing machine learning or analysis. 

### Some things to note
- My team continually gives me new data to play with, the data presented here is just a portion of what I am working with. 

- This is an ongoing project and I'm hoping to figure out how to automate and make my analysis more accurate this before I graduate so my team can hand it off to the next stats guy they see fit to work with.

-I will not be posting the data on here. I don't own the data. Contact me and I would be happy to talk about it with you.  

- Please message me if you have any suggestions, criticisms, questions, etc. I'm always open for discussion

![](thing.png)

### The Code and Data 
In this repository you'll find 3 different pieces of code. 
- Heatmap Script - I used a few plotting package
- Boxplot Script - The boxplots were all done on ggplot. I'm trying to figure how to automate this (not always possible with messy data)
- LASSO Script - This is my attemp (still a work in progress) to conduct variable selection and create a prediction model for the group I'm currently working for.
