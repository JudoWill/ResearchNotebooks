# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Learning Programming Through Examples

# <rawcell>

# This is the IPython notebook. You access it through the browser window however it is actually a program running on the super-computer on Will's desk. I have done some work to make sure that nothing you do here will screw up the computer. So feel free to experiment! This is a 'real' python interpreter and can do anything that python can do. Evaluate the next set of code-blocks (use 'shift-enter') to run the code in a cell.

# <codecell>

print 'Subraction', 500-300
print 'Multiplication', 500*300
print 'Raise to the Power', 500**30

# <rawcell>

# While having a $10,000 calculator is pretty cool, real programming involves using variables. Think of variables as just NAMES. Names can be anything WITHOUT SPACES, dashes, &, or other weird characters.

# <codecell>

this_is_a_number = 5
ThisIsAName = 'John Doe'
this_is_a_bigger_number = 500

# <rawcell>

# Then you can do anything with a variable that you can with a 'literal' (something you typed out).

# <codecell>

print 'Subtraction', this_is_a_number-this_is_a_bigger_number
print 'Multiplication', this_is_a_number*this_is_a_bigger_number

# <rawcell>

# You can also do interesting things. Like multiply a string (a variable with text) with a number. Can you guess what will happen?

# <codecell>

print ThisIsAName*this_is_a_number

# <rawcell>

# Oooooo ... concatenation! Think of it this way: 'What would you expect to happen if you had "5 times John Doe"'? What do you think would happen if you added a number and a string?

# <codecell>

print ThisIsAName+this_is_a_number

# <rawcell>

# Whoops,it doesn't understand. Python thinks we are trying to concatenate a string and a number. But numbers and strings aren't the same "Type" (duh). So it raises a TypeError and gives you a little description of what went wrong. This is called a Traceback ... when you're learning programming you'll see a lot of these. We'll go over them in more detail later ... but the important thing to know now is "Read them from the bottom to the top".

# <rawcell>

# Single varaibles are for chumps! List, Tuples (pronouced 'two-pulls'), sets, and dictionaries are ways to store multiple items in one name. Each has advantages and disadvantages.

# <codecell>

patient_ages =    (23,45,89,17,56,94) #Tuple, created using parentheses
patient_weights = [66,76,54,78,61,60] #List, created using brackets

# <rawcell>

# Items can be retrieved from a List or Tuple by thier index (position in the list). Just remember, computers start counting at ZERO!!!!

# <codecell>

print "Patient 4's age:", patient_ages[3]
print "Patient 4's weight:", patient_weights[3], "in kilos"
print 'In pounds:', patient_weights[3]*2.2

# <rawcell>

# There is one big difference between Lists and Tuples. Items in a list can be changed while items within a Tuple cannot be. This is called 'mutability'

# <codecell>

#This works
print 'Original:', patient_weights

patient_weights[4] = 5000
print 'After changing:', patient_weights

patient_weights[4] = 61
print 'After changing back', patient_weights


# <codecell>

#This does not
patient_ages[3] = 12
print patient_ages

# <rawcell>

# Dictionaries (pronounced and typed as 'dict') are UNORDERED mappings between keys and values. Just like a normal dictionary a key (kinda like a word) references one (and only one) value (kinda like a definition). Here are some examples.

# <codecell>

#This dictionary maps names to patient-numbers
names2index = {
                    'John':0,
                    'Jamie':1,
                    'Jack':2,
                    'Jill':3,
                    'Jess':4,
                    'Jeb':4,
}

print 'Jack is patient #', names2index['Jack']
print 'Jeb is patient #', names2index['Jeb']

#just remember, these are not 'in order' the way you normally think about dictionaries
print names2index
print sorted(names2index)

# <rawcell>

# So this lets us do some cool things:

# <codecell>

print "Jack's Age:", patient_ages[names2index['Jack']]
print "John's Weight in Kilos:", patient_weights[names2index['John']]
print "John's Weight in pounds:", patient_weights[names2index['John']]*2.2

# <rawcell>

# Lists, dicts and tuples allow us to do cool things called 'iteration'. This is one of the 'fundamental concepts' in programming, the idea of defining a set of rules for handeling an 'item' and then iterate through multiple 'items' doing the same thing. So as an example lets print the weight in kilos of each patient.

# <codecell>

for patient in names2index: #gets the name of each patient
    weight_kg = patient_weights[names2index[patient]]
    weight_lb = weight_kg*2.2
    print patient, 'weighs', weight_lb, 'in pounds'

# <rawcell>

# However keeping all of these things (Age, weight, name, etc) in seperate variables gets annoying. Especially as we start adding the hundreds of other things we could measure. So we're going to put them in a 'DataFrame'. Think of a DataFrame like an excel worksheet ... it has rows and columns and each position has a value (or is empty, something we'll cover later). However, 'python' doesn't have a native DataFrame (or anything like it); but thankfully people much smarter than anyone here have developed THOUSANDS of tools that can be accessed through the 'import' command. We'll cover them as we go; but in general there is a python module that can do almost anything you can think of!

# <codecell>

import antigravity

# <codecell>

#but seriously, the DataFrame object is a 'pandas' tool. You can see what else is there: http://pandas.pydata.org/
from pandas import DataFrame
# by importing in this way, you call it as DataFrame rather than pandas.DataFrame, so its a good way to save some space when writing code
# you can also do stuff like: import numpy as np (if you want to, but it can be considered bad style)

pat_data = DataFrame({
                      'Age':patient_ages,
                      'Weight':patient_weights
                        })
print pat_data

print 'Average age', pat_data['Age'].mean()
print 'Smallest weight', pat_data['Weight'].min()
print 'Oldest Age', pat_data['Age'].max()

# <codecell>

pat_data.T

# <rawcell>

# Once you 'grok' the concept of a DataFrame and all of the cool things you can do then you'll have almost limitless ability in dealing with 'tabular' data.

# <rawcell>

# For example, what if we wanted to find the number of patients with an Age > 50.

# <codecell>

print 'Make a "mask" of which patients have an Age > 50'

# a mask is basically a list of booleans (binary indexing)
old_mask = pat_data['Age']>50
print old_mask

print "Then we can 'sum' the True values"
print old_mask.sum()

# <rawcell>

# We can also use that 'mask' to index back into the DataFrame to answer a question like "Of all the patients over 50 what is thier average weight?"

# <codecell>

print "Old people's weight"
print pat_data['Weight'][old_mask]

# because python is mostly object-oriented programming, most things are methods rather than functions 
print "Old people's average weight in kg", pat_data['Weight'][old_mask].mean()
print "Old people's average weight in lbs", pat_data['Weight'][old_mask].mean()*2.2
print 

#A tildle (~) indicates 'not'
print "Young people's age"
print pat_data['Age'][~old_mask]

print "Young people's average weight in kg", pat_data['Weight'][~old_mask].mean()
print "Young people's average weight in lbs", pat_data['Weight'][~old_mask].mean()*2.2

# <rawcell>

# Hmm, it seems like young people weigh (on average) more then old people. Wouldn't it be cool if we could do a quick t-test to see if the difference is significant? As I said before, python has a module to do everything! Statistical functions like that are in the module "scipy.stats". You can see a full list of stat tests here: http://docs.scipy.org/doc/scipy/reference/stats.html (scroll about half-way down to find ttest_ind)

# <codecell>

from scipy.stats import ttest_ind

#this gives you a little pop-up showing the function definition (how to use it)
ttest_ind?

# <rawcell>

# There is a LOT of information in there. It takes practice to learn what the important things are but I'll walk you through this one.
# 
# "Definition: ttest_ind(a, b, axis=0, equal_var=True)"
# This tells us how to 'call' the function. It takes 2 'positional arguements' "a" and "b" and some optional arguements. We'll ignore them for now.
# 
# "Docstring:
# Calculates the T-test for the means of TWO INDEPENDENT samples of scores.
# 
# This is a two-sided test for the null hypothesis that 2 independent samples
# have identical average (expected) values. This test assumes that the
# populations have identical variances."
# 
# This tells us what the function does. All good libraries have some documentation telling what each function does.
# 
# "
# Parameters
# ----------
# "
# This tells what each input 'means'. So you can see the first two are 'array-like' which is just math/programming speak for "list of numbers". It tells us what the 'axis' arguement does, but it's just Greek now, so ignore it. It also tells us that we can tell the function whether we assume the populations have equal variance. For no we'll assume 'yes', which is the default.
# 
# "
# Returns
# -------
# t : float or array
#     The calculated t-statistic.
# prob : float or array
#     The two-tailed p-value.
# "
# This tells us that it returns two values, the t-statistic and the two-tailed p-value.
# 
# So we use it like this:

# <codecell>

old_weights = pat_data['Weight'][old_mask]
young_weights = pat_data['Weight'][~old_mask]
tstat, pval = ttest_ind(old_weights, young_weights)
print 'P-value:', pval

# <rawcell>

# Wow, pretty cool! We can say that (in this made-up data) old people have a statistically significantly (p<0.05) lower body weight then younger patients. But this is made up fake data. Would this be true in real data? Like the stuff we have from real patients? Well, lets see!

# <rawcell>

# I have done a lot of work to get data from redcap into python in a DataFrame so we can analyze questions like this. So you can load the redcap data into this python workspace and ask questions like this. How I did this is beyond the scope of this tutorial and there will be a little bit of 'magic' required to get it in. I'll explain as I go.

# <codecell>

from pandas import HDFStore #This is a file storage format for large collections of data

#This opens a particular version of the file for us to play with.
store = HDFStore('/home/will/HIVReportGen/Data/BaseRedcap/HIVAIDSGeneticAnalys_DATA_LABELS_2013-01-16_1211.hdf')

#This 'store' has multiple pieces of data within it
print store

#we want the redcap data (although the seq_data will be cool later
redcap_data = store['redcap']

#After we get the data out we 'close' the store, this is to make sure we can't accidentally save data BACK INTO the store.
#Think of it like 'taking out the jump-drive' once its copied to your desktop you don't need the jump-drive anymore.
store.close()

# <codecell>

print redcap_data

# <rawcell>

# From this you can see that there are 425 columns in the table, 115 of them are numbers (float), 171 are True/False (boolean) and 139 are something else (object) [usually these are strings]. But that still doesn't give us much information about WHAT is in there. But this will give us an idea:

# <codecell>

#print out the columns
print redcap_data.columns

# <rawcell>

# From that giant block of text we can see all the columns that we have. To answer out question we only need a small subset of them. So we'll pull out the Age and weight columns.

# <codecell>

real_pat_data = redcap_data[['Age', 'Weight']]
print real_pat_data

# <rawcell>

# Now you can start to see some of the hidden messyness of dealing with REAL data. You can see that not all of the data is there. We have 1419 entries only 1397 of them have an Age and only 544 of them have a Weight. But we'll work with what we have.

# <codecell>

real_old_mask = real_pat_data['Age']>50
print 'Old average weight:', real_pat_data['Weight'][real_old_mask].mean()
print 'Young average weight:', real_pat_data['Weight'][~real_old_mask].mean()

# <rawcell>

# Looks promissing. I wonder what the p-value is on that.

# <codecell>

tstat, pval = ttest_ind(real_pat_data['Weight'][real_old_mask], real_pat_data['Weight'][~real_old_mask])
print 'p-value:', pval

# <rawcell>

# Oooo, here we get into more of the problems with missing data. The t-test function doesn't know what to do with missing data (NaN) and so it choked and returned a NaN. So we need to think a second about what to do with missing data. It really depends on your situation, but in this one; the "right" answer is to just drop the missing data. And in this case we need to drop the patients who are missing an Age or Weight or BOTH. DataFrames have easy ways to deal with this: called "dropna". Lets look at the docstring.

# <codecell>

real_pat_data.dropna?

# <rawcell>

# Again, lots of stuff but only a little is important.
# 
# "Definition: real_pat_data.dropna(self, axis=0, how='any', thresh=None, subset=None)"
# How to call the function. The 'self' arguement is automatically added, due to "magic" you don't need to add it, as you get more programming experience it will make more sense as to why, but for right now, just know "skip anything that says 'self'"
# 
# "
# Docstring:
# Return object with labels on given axis omitted where alternately any
# or all of the data are missing"
# This tells us that it will drop items from the dataset if there are missing values.
# 
# 
# In the paramenters we need to pay attention:
# 
# axis : {0, 1}, or tuple/list thereof
#     Pass tuple or list to drop on multiple axes
# This lets us change which 'axis' it will work on. Think of axis Zero as the 'rows' and axis One as the 'columns'. In our case we want to drop ROWS which have missing values not COLUMNS. So we'll use "axis = 0".
# 
# how : {'any', 'all'}
#     any : if any NA values are present, drop that label
#     all : if all values are NA, drop that label
# This lets us control how stringent the criteria is. Will it drop a row if ANY are missing or if ALL are missing.
# 
# thresh : int, default None
#     int value : require that many non-NA values
# If ANY or ALL don't work we could use an integer number of missing values. For example drop if more than 5 columns are empty in a 50 column DataFrame.
# 
# subset : array-like
#     Labels along other axis to consider, e.g. if you are dropping rows
#     these would be a list of columns to include
# This is for when we only want it to consider a SUBSET of the columns for checking missing values.
# 
# So in our case we want to drop the rows of the dataset inwhich any of the columns are empty.

# <codecell>

full_data = real_pat_data.dropna(axis = 0, how = 'any')
print full_data

# <rawcell>

# So, we have 543 rows in which we have both an Age and a Weight, now we can do something.

# <codecell>

tstat, pval = ttest_ind(full_data['Weight'][real_old_mask], full_data['Weight'][~real_old_mask])
print 'p-value:', pval

# <rawcell>

# Wow, it worked! Old people do weigh less the young people (p<0.05). We also hiddenly did something else REALLY cool, that you didn't even see, did you?
# 
# 
# We used the 'mask' (real_old_mask) calculated on the original dataset (real_pat_data) on the filtered dataset (full_data). If you look at the "Index" in the print-out from each dataset you can get a hint as to how this happened.

# <codecell>

#look at the first 10 items in the index (think row-labels)
print 'Large Data row labels:', real_old_mask.index[:10]
print 'Smaller Data row labels:', full_data.index[:10]

# <rawcell>

# Pandas was smart enough to go through and match the rows which have the same row labels. This can be good and bad. It makes life much easier but can sometimes lead to errors if you're not careful!

# <rawcell>

# But there's an error in our statistics logic ... the code is right, but our analysis is wrong. You know that we have multiple patients dataset who have multiple visits (that whole longitudinal thingy). So some patients are counted more then once and thier ages/weights are likely to be similar. If our cohort has a bias (intentional or accidental) in the age or weight selection for returning patients then it will bias our test. So we should probably consider only one Age and one Weight for each patient. There are COUNTLESS ways to do this, I'm going to show a few and use them as ways to show you other cool Python-y things.

# <rawcell>

# First lets look at each patient and find thier average age and average weight. To do that we're going to use the 'groupby' function to group rows by Patient ID and then calculate the average age and weight.

# <codecell>

#pull out the data
data = redcap_data[['Patient ID', 'Age', 'Weight']]
print 'raw data\n'
print data

# groupby is by far the best function so far - this will save SO MUCH TIME
mean_data = data.groupby('Patient ID').mean()
print '\nMean-ed data\n'
print mean_data

# Drop the NaN values
print '\nFull without missing values\n'
full_mean_data = mean_data.dropna(axis = 0, how = 'any')
print full_mean_data

# <rawcell>

# You can see that it also changed the 'index' of the mean_data DataFrame. It is now indexed by Patient ID "Index: 507 entries, A0001 to A0514" and we can also see that many of the patients do not have Ages and Weights for ANY of thier visits. So we'll work with the data we got and see what happens.

# <codecell>

tstat, pval = ttest_ind(full_mean_data['Weight'][real_old_mask], full_mean_data['Weight'][~real_old_mask])
print 'p-value:', pval

# <rawcell>

# Whoops, an error! Hmmm, what's happening? It says: "ValueError: cannot index with vector containing NA / NaN values". That doesn't make much sense ... we dropped all the NaN values from "full_mean_data" and we've used this "real_old_mask" variable before and it worked just fine. What could it mean? Let's deconstruct the code above and see where the error is.

# <codecell>

old_weights = full_mean_data['Weight'][real_old_mask]
young_weights = full_mean_data['Weight'][~real_old_mask]
tstat, pval = ttest_ind(old_weights, young_weights)
print 'p-value:', pval

# <rawcell>

# Still have an error (duh, I wrote the same code, just more 'verbose'). But if we read it bottom-to-top then you can see something. The error seems to be when we try to pull out the weights. Hmmmm. Let's look at the data and see if there's something fishy.

# <codecell>

print 'Data'
print full_mean_data.head() #prints only the first 5 rows, useful for checking sanity
# can also do .tail() for the last 5

print
print 'Mask'
print real_old_mask.head()

# <rawcell>

# Eureaka! Do you see it? The first column in both tables is the Index. The Data is indexed by "Patient ID" and the Mask is indexed by numbers. And those numbers refer to the rows in the 'real_pat_data' matrix ... so they can't be mapped to Patient ID's. So we just need to make a new mask! Remember how I said above that Pandas is smart ... well sometimes its just smrt.

# <codecell>

pat_old_mask = full_mean_data['Age']>50
print pat_old_mask.head()

# <rawcell>

# There we go! Now we know which patients have an average age > 50! Now we can go back to our original question.

# <codecell>

tstat, pval = ttest_ind(full_mean_data['Weight'][pat_old_mask], full_mean_data['Weight'][~pat_old_mask])
print 'p-value:', pval

# <rawcell>

# Interesting ... when we correct for the fact that we have patients with multiple visits (and therefore the age and weights from each row are not TRULY independent) we lose significance. Maybe the 'average' weight and 'average' age is not the correct way to think about it. It would make more sense to look at a single visit for each patient and take the Age and Weight from that visit. But which Visit? To be 'fair' we need to do it one the same visit for each patient. One could make a case for any particular visit, so I'm going to say thier 'initial visit' (ROO) since everyone has an ROO.
# 
# A little note about the database. Its niave to believe that its 'sorted' in any particular order. So make sure any of your tests don't rely on that, always check first.

# <codecell>

#so lets pull out the visit data along with Patient ID, Age, Weight.
data = redcap_data[['Patient ID', 'Patient visit number', 'Age', 'Weight']]
print data
print data.head()

# <rawcell>

# Notice something? Not every visit has a Visit Number.
# There are many ways to do this, but I'll show the simplest.

# <codecell>

first_visit_mask = data['Patient visit number']=='R00'
first_data = data[first_visit_mask].dropna(axis = 0, how = 'any')
print first_data.head()

# <rawcell>

# Looks like it worked. Let's see what the result is.

# <codecell>

first_old_mask = first_data['Age'] > 50
tstat, pval = ttest_ind(first_data['Weight'][first_old_mask], first_data['Weight'][~first_old_mask])
print 'p-value:', pval

# <rawcell>

# Wow! Not even CLOSE! That seems weird doesn't it? When I look at the average weight across all visits and compare young/old I'm almost significant but when I look only at the intake visit its way off. I wonder what's happening here. Let's look at the data. We can use some plotting to see what's happening. Why don't we look at the weight at each visit. We'll use a boxplot to do this.

# <codecell>

data.boxplot?
#you should be able to read these now!

# <codecell>

data.boxplot(column = 'Age', by = 'Patient visit number')
data.boxplot(column = 'Weight', by = 'Patient visit number')

# <rawcell>

# It seems like the R09 and R10 visits have different variances and much lower medians. Maybe if we exclude those we'll get answers that match better.

# <codecell>

rm_mask = (data['Patient visit number'] == 'R09') | (data['Patient visit number'] == 'R10')
early_data = data[~rm_mask].dropna(axis = 0, how = 'any')
print early_data.head()

grouped_early_data = early_data.groupby('Patient ID').mean()
early_mask = grouped_early_data['Age'] > 50
tstat, pval = ttest_ind(grouped_early_data['Weight'][early_mask], grouped_early_data['Weight'][~early_mask])
print 'p-value:', pval

# <rawcell>

# Hmm, That didn't seem to help. It still looks like the average weight of patients less than 50 yrs-old is different from the average weight of patients older than 50. But why is is it SO different from a result where we calculate based on the intake visit? Lets look again.

# <rawcell>

# Hmmmm, I don't see anything else useful here. The variance of the age/weight decreases as we look at higher visits but that's to be expected (less patients make it to the later visits). Maybe its something about the 'averaging'. Maybe young people and old people have different variances in thier weights. We could look at that by calculating the Standard Deviation for the weight of each patient.

# <codecell>

pat_weight_std = data.groupby('Patient ID', as_index = False).aggregate({'Weight':'std'})
#I'll explain as_index later.
print pat_weight_std.describe()

# <rawcell>

# Hmm, so some people have a HUGE variation in weight (67 lbs std is a lot). Lets see if that is correlated with Age. So lets put it in with the Age at thier first visit.

# <codecell>

first_data['Age-std'] = pat_weight_std
print first_data

# <rawcell>

# Hmmm, it won't let us put them together saying "Length of values does not match length of index" ... Implying that the two items don't match properly. Since this is a common issue, Pandas has a way to deal with this. Called, "merge".

# <codecell>

from pandas import merge
merge?

# <codecell>

mdata = merge(first_data, pat_weight_std,
                left_on = 'Patient ID', right_on = 'Patient ID')
print mdata.head()

# <rawcell>

# Hmmm, what's going on with that Weight_x, Weight_y? It comes when doing a 'merge' and the two DataFrames have columns with the same name.

# <codecell>

print first_data.columns
print pat_weight_std.columns
print 'See, Weight is in both of them.'

# <rawcell>

# So we can rename the column in one and then it will work fine.

# <codecell>

pat_weight_std = pat_weight_std.rename(columns={'Weight':'Weight-std'})

mdata = merge(first_data, pat_weight_std,
                left_on = 'Patient ID', right_on = 'Patient ID')
print mdata.head()

# <rawcell>

# Now we can go back to our earlier question? Does the weight vary more or less with age and does that account for our earlier observations.

# <codecell>

mdata[['Age', 'Weight-std']].scatter()

# <rawcell>

# Hmm, there's no scatter plot in DataFrame. I wonder if there's another package that does it? Why don't you google 'scatter plot in python'. What is the package in the first few results?

# <codecell>

from matplotlib import pyplot

# <codecell>

pyplot.scatter(mdata['Age'], mdata['Weight-std'])
pyplot.xlabel('Age')
pyplot.ylabel('Weight-std')

# <rawcell>

# Hmmm, I don't see much of a trend. However, that one patient is WAY out there in terms of weight-control and they're right on the cusp of 50. What happens if we drop them from the analysis. First we have to find them.

# <codecell>

print mdata[mdata['Weight-std'] > 60]

# <rawcell>

# If we wanted to drop a patient from a dataset what do you think the method will be?

# <codecell>

full_mean_data.drop?

# <codecell>

tstat, pval = ttest_ind(full_mean_data['Weight'][pat_old_mask], full_mean_data['Weight'][~pat_old_mask])
print 'Original p-value:', pval

trimmed_data = full_mean_data.drop(['A0502'], axis = 0)
tstat, pval = ttest_ind(trimmed_data['Weight'][pat_old_mask], trimmed_data['Weight'][~pat_old_mask])
print 'Trimmed p-value:', pval

# <rawcell>

# Hmm, that seemed to help fix the p-value, but not much. I wonder if there's something else we can find.

# <codecell>


