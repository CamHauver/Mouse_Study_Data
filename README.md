# Mouse_Study_Data
Pymaceuticals, Inc. mouse drug regimen study data analysis; data visualizations and summary statistics for technical report of clinical study

#Original code sources: #Data Analytics course pymaceuticals_starter code from Module 5 Challenge files #Data Analytics course instructor, Andrew Hoang's, speed run Zoom recording for Module 5 Challenge

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# - See statistics summary (summary_statistics_df) and various data visualizations below for visual summary of various aspects of this study's results
#  

# In[1]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

#print(mouse_metadata.head())
#print(mouse_metadata.shape)
#print(study_results.head())
#print(study_results.shape)

# Combine the data into a single DataFrame
mouse_study_data = pd.merge(
    study_results, mouse_metadata, how="left", on=["Mouse ID"]
)
#print(mouse_study_data.shape)

# Display the data table for preview
mouse_study_data.head()


# In[2]:


# Checking the number of mice.
len(mouse_study_data["Mouse ID"].unique())


# In[3]:


# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint.

dup_mouse_id = mouse_study_data[mouse_study_data.duplicated(subset=["Mouse ID", "Timepoint"])]["Mouse ID"].unique()
dup_mouse_id


# In[4]:


# Optional: Get all the data for the duplicate mouse ID.
dup_mouse_id_all_data = mouse_study_data[mouse_study_data["Mouse ID"] == "g989"]
dup_mouse_id_all_data


# In[5]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_mouse_study_data = mouse_study_data[mouse_study_data["Mouse ID"].isin(dup_mouse_id) == False]
clean_mouse_study_data


# In[6]:


# Checking the number of mice in the clean DataFrame.
len(clean_mouse_study_data["Mouse ID"].unique())


# ## Summary Statistics

# In[7]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen:
# mean, median, variance, standard deviation, and SEM of the tumor volume.

means = clean_mouse_study_data.groupby("Drug Regimen").mean(numeric_only=True)["Tumor Volume (mm3)"]
medians = clean_mouse_study_data.groupby("Drug Regimen").median(numeric_only=True)["Tumor Volume (mm3)"]
variances = clean_mouse_study_data.groupby("Drug Regimen").var(numeric_only=True)["Tumor Volume (mm3)"]
stdevs = clean_mouse_study_data.groupby("Drug Regimen").std(numeric_only=True)["Tumor Volume (mm3)"]
sems = clean_mouse_study_data.groupby("Drug Regimen").sem(numeric_only=True)["Tumor Volume (mm3)"]

# Assemble the resulting series into a single summary DataFrame.
summary_statistics_df = pd.DataFrame({
    "Mean Tumor Volume (mm3)": means,
    "Median Tumor Volume (mm3)": medians,
    "Tumor Volume (mm3) Variance": variances,
    "Tumor Volume (mm3) Standard Deviation": stdevs,
    "Tumor Volume (mm3) SEM": sems,
})
summary_statistics_df


# In[8]:


# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line
clean_mouse_study_data.groupby("Drug Regimen").agg({"Tumor Volume (mm3)":["mean", "median", "var", "std", "sem"]})


# ## Bar and Pie Charts

# In[9]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
mouse_id_count = clean_mouse_study_data["Drug Regimen"].value_counts()

mouse_id_count.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.ylabel("# of Observed Mouse Timepoints")


# In[10]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
plt.bar(mouse_id_count.index.values, mouse_id_count.values)
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("# of Observed Mouse Timepoints")
plt.show


# In[19]:


# Generate a pie chart, using Pandas, showing the distribution of unique female versus male mice used in the study

# Get the unique mice with their gender
mouse_sex_counts = clean_mouse_study_data.Sex.value_counts()

# Make the pie chart
mouse_sex_counts.plot(kind="pie", autopct="%1.1f%%")
plt.ylabel("Sex Proportion")


# In[20]:


# Generate a pie chart, using pyplot, showing the distribution of unique female versus male mice used in the study

# Get the unique mice with their gender
mouse_sex_counts = clean_mouse_study_data.Sex.value_counts()

# Make the pie chart
plt.pie(mouse_sex_counts.values, labels = mouse_sex_counts.index.values, autopct="%1.1f%%")
plt.ylabel("Sex Proportion")


# ## Quartiles, Outliers and Boxplots

# In[13]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
max_tumor_size_data = clean_mouse_study_data.groupby(["Mouse ID"])["Timepoint"].max()
max_tumor_size_data = max_tumor_size_data.reset_index()

# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
mouse_study_data_final = max_tumor_size_data.merge(clean_mouse_study_data, on=["Mouse ID", "Timepoint"], how="left")
mouse_study_data_final


# In[14]:


# Put treatments into a list for for loop (and later for plot labels)
treatment_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_data = []

# Calculate the IQR and quantitatively determine if there are any potential outliers.
for drug in treatment_list:

    # Locate the rows which contain mice on each drug and get the tumor volumes
    mouse_study_data_final_tumor_vol = mouse_study_data_final.loc[mouse_study_data_final["Drug Regimen"] == drug, "Tumor Volume (mm3)"]

    # add subset
    tumor_vol_data.append(mouse_study_data_final_tumor_vol)

    # Determine outliers using upper and lower bounds
    quartiles = mouse_study_data_final_tumor_vol.quantile([0.25, 0.5, 0.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq - lowerq
    lower_bound = lowerq - (1.5 * iqr)
    upper_bound = upperq + (1.4 * iqr)

    outliers = mouse_study_data_final_tumor_vol.loc[(mouse_study_data_final_tumor_vol < lower_bound) | (mouse_study_data_final_tumor_vol > upper_bound)]
    print(f"{drug}'s potential outliers {outliers}")


# In[15]:


# Generate a box plot that shows the distribution of the tumor volume for each treatment group.
orange_out = dict(markerfacecolor="red", markersize=15)
plt.boxplot(tumor_vol_data, labels = treatment_list, flierprops=orange_out)
plt.ylabel("Final Tumor Volume (mm3")
plt.show()


# ## Line and Scatter Plots

# In[16]:


# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
capomulin_data = clean_mouse_study_data[clean_mouse_study_data["Drug Regimen"] == "Capomulin"]
mouse_l509_capomulin = capomulin_data[capomulin_data["Mouse ID"] == "l509"]
plt.plot(mouse_l509_capomulin["Timepoint"], mouse_l509_capomulin["Tumor Volume (mm3)"])
plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor Volume (mm3)")
plt.title("Capomulin treatment of mouse l509")


# In[21]:


# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen
capomulin_data_mean = capomulin_data.groupby(["Mouse ID"]).mean(numeric_only=True)
capomulin_data_mean

plt.scatter(capomulin_data_mean["Weight (g)"], capomulin_data_mean["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.title("Capomulin Treatment Average Tumor Volume vs Mouse Weight")


# ## Correlation and Regression

# In[22]:


# Calculate the correlation coefficient and a linear regression model
# for mouse weight and average observed tumor volume for the entire Capomulin regimen
correlation = st.pearsonr(capomulin_data_mean["Weight (g)"], capomulin_data_mean["Tumor Volume (mm3)"])
print(f"The correlation between mouse weight and average tumor volume is {round(correlation[0],2)}")

model = st.linregress(capomulin_data_mean["Weight (g)"], capomulin_data_mean["Tumor Volume (mm3)"])
slope = model[0]
intercept = model[1]
y_values = capomulin_data_mean["Weight (g)"] * slope + intercept
plt.scatter(capomulin_data_mean["Weight (g)"], capomulin_data_mean["Tumor Volume (mm3)"])
plt.plot(capomulin_data_mean["Weight (g)"], y_values, color="red")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.title("Capomulin Treatment Average Tumor Volume vs Mouse Weight")
plt.show()
