# pegRNA_rank
Deep learning MLP model to rank pegRNA sequence

The code is structured to take input from a user for two arguments (First: Target Region Amplicon, Second: Desired Insert). Then the program will output a csv file containing the predicted score for each sequence named "Result.csv" within the same directory

# Setup Instructions
Prerequisites: Make sure you have Python>=3.6 and PyTorch>=1.5 installed. Then, install dependencies with:
Then run
pip install -r requirements.txt

Next, either clone and add the base directory of the repository to your PYTHONPATH

Then run
python rank.py

