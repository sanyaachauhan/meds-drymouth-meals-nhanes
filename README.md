# Meds, Dry Mouth, and Meals: What NHANES Data Tells Us About Older Adults

**Author:** Sanya Chauhan | MSc Health Sciences, University of Lucerne  
**Background:** Dentist (BDS)  

## What this is about
In the dental clinic, one of the most common things observed in older adults isn't just tooth decay, it's severe dry mouth (xerostomia). Patients come in taking multiple daily medications for blood pressure, allergies, sleep and they talk about how chewing feels like swallowing cardboard, food loses its taste and eating becomes a chore. 

I wanted to step outside the dental chair and look at the bigger picture: Does medication induced dry mouth actually mess up a person's nutrition at a national scale? And is there a specific tipping point where things fall apart?

To find out, I built an R pipeline using public health data from the US government to map out how medication burden connects to oral health and daily caloric intake.

---

## What data and tools I used
* **Dataset:** CDC NHANES (National Health and Nutrition Examination Survey) spanning 2011–2018. I combined 4 continuous survey cycles to get a solid sample size of roughly **39,000 participants**.
* **Medication Scoring:** Wrote a custom script to scan prescription drug names and calculate an **Anticholinergic Cognitive Burden (ACB)** score giving points based on how badly a drug dries out the mouth.
* **Survey Weights:** Applied official CDC complex survey weights so the numbers actually reflect the broader US population rather than a biased sample.
* **Analysis & Machine Learning:** Ran standard survey regressions (GLMs) and a survey weighted **XGBoost** model with **SHAP** plots to spot non linear trends.

---

## What I found
1. **The Tipping Point:** Medication side effects don't just scale in a straight line. Once a patient's cumulative ACB score hits **ACB ≥ 2** (meaning multiple drying drugs), their daily caloric intake takes a sharp, compounding plunge.
2. **The "Double Hit" Group:** About **0.87% of US adults (~2.47 million people)** are trapped in a high risk overlap: they have both a high dry mouth drug burden *and* poor oral health at the same time.
3. **Why it matters:** Instead of triggering false alarms for every single person with dry mouth, targeting this 0.87% high risk cohort could give doctors and pharmacists a practical way to step in before severe nutritional frailty sets in.

---

## Code & Files in this Repo
* `01_data_extraction.R`: The main analysis script (pulls data, cleans columns, builds the ACB score, runs models and SHAP plots).
* `Population_Risk_Exposure_Summary.csv`: The exported population risk breakdown.
* `Figure1_Global_Feature_Hierarchy.png`, `Figure2_Cascade_Interaction.png`, `Figure3_Tipping_Point.png`: Visual outputs showing feature importance, interaction curves, and the caloric tipping point.

*Built using R (`tidyverse`, `nhanesA`, `survey`, `xgboost`, `shapviz`).*
