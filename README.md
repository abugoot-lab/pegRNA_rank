# pegRNA_rank
Deep learning MLP model to rank pegRNA sequence

User Instruction:
1. When prompted by the program, the user has to input a target region of interest (amplicon) that they would like to perform insertion first. 
2. Then the program will ask whether user would like to specificy a custom insertion sequence, that the user can respond by "Y" or "N" (representing yes or no respectively). 
3. If the user select yes, then the user will be asked to input a customized insertion sequence.
4. Otherwise, the program will automatically run by using the default insertion sequences in the database
5. Lastly, the program will output a csv file containing the predicted score for each atgRNA sequence named "Result.csv" within the same directory

# Setup Instructions
Prerequisites: Make sure you have Python>=3.6 and PyTorch>=1.5 installed. Then, install dependencies with:
```
pip install -r requirements.txt

```

Next, clone and add the base directory of the repository to your PYTHONPATH

Then run:
```
python rank.py
```
