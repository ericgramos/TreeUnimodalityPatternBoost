from tokenizers import Tokenizer
import os
import argparse

# Set up argument parser for command line options
parser = argparse.ArgumentParser()
parser.add_argument('--input-file', '-i', type=str, default='out/out.txt', help="input file with tokenized data")
parser.add_argument('--output-file', '-o', type=str, default='transformer-output-decoded.txt', help="output file for decoded text")
parser.add_argument('--tokenizer-dir', '-t', type=str, default='tokenizer_data', help="directory where tokenizer is saved")
args = parser.parse_args()

# Load the tokenizer from the saved file
tokenizer_path = os.path.join(args.tokenizer_dir, "tokenizer.json")
if not os.path.exists(tokenizer_path):
    raise FileNotFoundError(f"No tokenizer found at {tokenizer_path}. Please check the path and try again.")

tokenizer = Tokenizer.from_file(tokenizer_path)

# Function to decode a single line of token ids
def decode_tokens(token_line):
    # Remove the 'V' prefix and convert to integers
    #print(token_line)
    token_ids = [int(token[1:]) for token in token_line.split(',')]
    # Decode the token ids to text
    
    return tokenizer.decode(token_ids).replace(" ","")

# Process the input file
if os.path.exists(args.input_file):
    with open(args.input_file, 'r') as file:
        tokenized_lines = file.readlines()

    # Decode each line and collect the results
    decoded_text = [decode_tokens(line.strip()) for line in tokenized_lines if len(line) > 1]

    # Write the decoded text to the output file
    with open(args.output_file, 'w') as file:
        for line in decoded_text:
            file.write(line + '\n')

    print(f"Decoding complete. Check the output in {args.output_file}")
else:
    print(f"Error: The file {args.input_file} does not exist.")

