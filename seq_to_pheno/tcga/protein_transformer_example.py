import torch
import transformers
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from transformers import BertTokenizer, BertModel
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
import numpy as np
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Load the dataset with error handling
try:
    df = pd.read_csv('seq_to_pheno/tcga/data/protein_sequences_metadata.tsv', sep='\t')
    logger.info("Dataset loaded successfully.")
except FileNotFoundError as e:
    logger.error(f"File not found: {e}")
    raise
except pd.errors.EmptyDataError as e:
    logger.error(f"Empty data error: {e}")
    raise
except Exception as e:
    logger.error(f"An error occurred while loading the dataset: {e}")
    raise

# Extract relevant columns
try:
    sequences = df['mutated_protein'].values
    transcript_ids = df['transcript_id'].values
    survival_times = df['Donor Survival Time'].values
    logger.info("Relevant columns extracted successfully.")
except KeyError as e:
    logger.error(f"Key error: {e}")
    raise

# Split the dataset into training and validation sets
try:
    train_sequences, val_sequences, train_labels, val_labels, train_ids, val_ids = train_test_split(
        sequences, survival_times, transcript_ids, test_size=0.2, random_state=42
    )
    logger.info("Dataset split into training and validation sets.")
except Exception as e:
    logger.error(f"An error occurred during dataset splitting: {e}")
    raise

# Define custom dataset class
class ProteinDataset(Dataset):
    def __init__(self, sequences, labels, tokenizer, max_length):
        self.sequences = sequences
        self.labels = labels
        self.tokenizer = tokenizer
        self.max_length = max_length

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        label = self.labels[idx]
        
        try:
            # Tokenize the sequence
            encoding = self.tokenizer(
                sequence,
                add_special_tokens=True,
                max_length=self.max_length,
                padding='max_length',
                truncation=True,
                return_tensors='pt'
            )
        except Exception as e:
            logger.error(f"An error occurred during tokenization: {e}")
            raise

        input_ids = encoding['input_ids'].squeeze()
        attention_mask = encoding['attention_mask'].squeeze()

        return {
            'input_ids': input_ids,
            'attention_mask': attention_mask,
            'labels': torch.tensor(label, dtype=torch.float)
        }

# Initialize the tokenizer
try:
    tokenizer = BertTokenizer.from_pretrained('Rostlab/prot_bert')
    logger.info("Tokenizer initialized.")
except Exception as e:
    logger.error(f"An error occurred while initializing the tokenizer: {e}")
    raise

# Define hyperparameters
MAX_LENGTH = 256
BATCH_SIZE = 2
LEARNING_RATE = 1e-5
EPOCHS = 5

# Create datasets and dataloaders
train_dataset = ProteinDataset(train_sequences, train_labels, tokenizer, MAX_LENGTH)
val_dataset = ProteinDataset(val_sequences, val_labels, tokenizer, MAX_LENGTH)

train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False)

# Define the model
class ProteinBERTRegressor(nn.Module):
    def __init__(self):
        super(ProteinBERTRegressor, self).__init__()
        try:
            self.bert = BertModel.from_pretrained('Rostlab/prot_bert')
            self.regressor = nn.Linear(self.bert.config.hidden_size, 1)
            logger.info("Model initialized.")
        except Exception as e:
            logger.error(f"An error occurred while initializing the model: {e}")
            raise

    def forward(self, input_ids, attention_mask):
        try:
            outputs = self.bert(input_ids=input_ids, attention_mask=attention_mask)
            pooler_output = outputs.pooler_output
            output = self.regressor(pooler_output)
            return output
        except Exception as e:
            logger.error(f"An error occurred during forward pass: {e}")
            raise

# Initialize the model, loss function, and optimizer
device = torch.device('cuda' if torch.cuda.is_available() else 'mps' if torch.backends.mps.is_available() else 'cpu')
model = ProteinBERTRegressor().to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)

# Training loop
train_losses = []
gradient_accumulation_steps = 4
val_losses = []
for epoch in range(EPOCHS):
    model.train()
    total_loss = 0
    for batch_idx, batch in enumerate(train_loader):
        input_ids = batch['input_ids'].to(device)
        attention_mask = batch['attention_mask'].to(device)
        labels = batch['labels'].to(device)

        optimizer.zero_grad()
        try:
            outputs = model(input_ids, attention_mask)
            loss = criterion(outputs.squeeze(), labels)
            loss = loss / gradient_accumulation_steps
            loss.backward()

            if (batch_idx + 1) % gradient_accumulation_steps == 0:
                optimizer.step()
                optimizer.zero_grad()
        except Exception as e:
            logger.error(f"An error occurred during training: {e}")
            raise

        total_loss += loss.item()

    avg_train_loss = total_loss / len(train_loader)
    train_losses.append(avg_train_loss)
    logger.info(f"Epoch {epoch+1}/{EPOCHS}, Training Loss: {avg_train_loss:.4f}")

    # Validation loop
    model.eval()
    val_loss = 0
    val_predictions = []
    val_labels_list = []
    val_ids_list = []
    with torch.no_grad():
        for batch_idx, batch in enumerate(val_loader):
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            labels = batch['labels'].to(device)

            try:
                outputs = model(input_ids, attention_mask)
                loss = criterion(outputs.squeeze(), labels)
                val_loss += loss.item()

                val_predictions.extend(outputs.squeeze().tolist())
                val_labels_list.extend(labels.tolist())
                val_ids_list.extend(val_ids[batch_idx * BATCH_SIZE:(batch_idx + 1) * BATCH_SIZE])
            except Exception as e:
                logger.error(f"An error occurred during validation: {e}")
                raise

    avg_val_loss = val_loss / len(val_loader)
    val_losses.append(avg_val_loss)
    logger.info(f"Epoch {epoch+1}/{EPOCHS}, Validation Loss: {avg_val_loss:.4f}")

# Plot training and validation loss
plt.plot(range(1, EPOCHS + 1), train_losses, label='Training Loss')
plt.plot(range(1, EPOCHS + 1), val_losses, label='Validation Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.title('Training and Validation Loss')
plt.legend()
plt.show()

# Identify transcripts with highest and lowest relationship to Donor Survival Time
sorted_indices = np.argsort(val_predictions)
lowest_transcripts = [val_ids_list[i] for i in sorted_indices[:5]]
highest_transcripts = [val_ids_list[i] for i in sorted_indices[-5:]]

logger.info("Transcripts with lowest predicted survival times:")
for transcript in lowest_transcripts:
    logger.info(transcript)

logger.info("\nTranscripts with highest predicted survival times:")
for transcript in highest_transcripts:
    logger.info(transcript)

# Explanation for implementing the model after training
# To implement the model after training, you can save the model state and use it for inference as follows:

# Save the model
try:
    model_path = 'seq_to_pheno/tcga/data/protein_bert_regressor.pth'
    torch.save(model.state_dict(), model_path)
    logger.info("Model saved successfully.")
except Exception as e:
    logger.error(f"An error occurred while saving the model: {e}")
    raise

# Load the model for inference
def load_model(model_path):
    model = ProteinBERTRegressor()
    try:
        model.load_state_dict(torch.load(model_path, map_location=device))
        model.to(device)
        model.eval()
        logger.info("Model loaded successfully.")
    except Exception as e:
        logger.error(f"An error occurred while loading the model: {e}")
        raise
    return model

# Example of using the model for inference
def predict(model, sequence):
    try:
        encoding = tokenizer(
            sequence,
            add_special_tokens=True,
            max_length=MAX_LENGTH,
            padding='max_length',
            truncation=True,
            return_tensors='pt'
        )

        input_ids = encoding['input_ids'].to(device)
        attention_mask = encoding['attention_mask'].to(device)

        with torch.no_grad():
            output = model(input_ids, attention_mask)
        return output.item()
    except Exception as e:
        logger.error(f"An error occurred during prediction: {e}")
        raise

# Load the model and make a prediction
try:
    model = load_model(model_path)
    example_sequence = "MDDYMVLRMIGEGSFGRALLVQHESSNQMFAMKEIRLPKSFSNTQN"
    predicted_survival_time = predict(model, example_sequence)
    logger.info(f"Predicted Survival Time: {predicted_survival_time:.2f}")
except Exception as e:
    logger.error(f"An error occurred during model prediction: {e}")