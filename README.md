# Stage Dependent Modeling of Cervical Cancer from Gene Expression Data (Machine Learning on Biological Data using R)

## Project Description:
I'm being funded to help the [Wallace Lab](https://wallacelabksu.weebly.com/) at Kansas State University. Specifically, I'm applying [LASSO](https://www.jstor.org/stable/2346178?seq=1), which is an efficient machine learning model used for feature selection and regularization. 

## Project Inspiration
Prophylactic vaccines can prevent the human papillomavirus (HPV) infections that cause cervical cancers (CaCx). Yet, CaCx still kills a woman every 90 seconds. Vaccine hesitancy hinders the herd immunity and broad protection in wealthy countries. Some of the hesitancy comes from the very low rates of CaCx thanks to widespread frequent pap smear screening that are quite effective at providing early interventions that prevent these malignancies. Indeed, the 14 million new HPV infections each year in the US, only result in 13,170 cases of CaCx. In developing countries, where the majority of CaCx deaths occur, fiscal and infrastructural constraints further any hope for a successful vaccination campaign and present significant barriers to pap smear screening. 
While these screening modalities are quite effective at detecting pre-malignant lesions, they could be done less frequently if we could identify women who are at higher risk of progression. Most HPV infections and many premalignant lesions would be cleared by the host immune response without ever posing a cancer risk. However, since there are no established biomarkers for triaging women based on the relative likelihood of progression, paps smears have become routine.  If women at greater risk could be identified, paps smear screening could be optimized by less frequent screening of low risk cohorts. In higher socioeconomic groups, this would decrease the stress and physical discomfort associated with the exam. The potential rewards in low and middle-income nations would be far greater. In these communities, reducing the number of interventions would allow limited resources to be better allocated to people with the greatest need.

### LASSO
I'm using Multinomial LASSO which is found in the GLMNET package description [here](https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html). It is specifically under the Multinomial Models section of the webpage. 

### Why LASSO?
LASSO is great for parsimonious modeling when dealing with high-dimesional n<<p datasets because of the shrinkage it conducts (feature selection). This helps us deal with overfitting. 

## Workflow
![]()
