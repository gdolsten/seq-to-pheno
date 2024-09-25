import os
from datetime import datetime
from huggingface_hub import HfApi, Repository
from typing import Optional
import os
import yaml
from typing import Dict, Any

class DatasetCard:
    def __init__(self, template_path: str):
        with open(template_path, 'r') as f:
            self.template = f.read()

    def generate(self, data: Dict[str, Any]) -> str:
        """Generate the dataset card content by filling in the template with provided data."""
        return self.template.format(**data)

    def save(self, output_path: str, data: Dict[str, Any]):
        """Save the generated dataset card to the specified file."""
        content = self.generate(data)
        with open(output_path, 'w') as f:
            f.write(content)

    def write_to_string(self, data: Dict[str, Any]) -> str:
        """Generate the dataset card and return it as a string."""
        return self.generate(data)

    def prepend(self, prepend_text: str, data: Dict[str, Any]) -> str:
        """Prepend text to the generated dataset card."""
        content = self.generate(data)
        return prepend_text + '\n' + content

    def postpend(self, postpend_text: str, data: Dict[str, Any]) -> str:
        """Append (postpend) text to the generated dataset card."""
        content = self.generate(data)
        return content + '\n' + postpend_text

    def prepend_and_postpend(self, prepend_text: str, postpend_text: str, data: Dict[str, Any]) -> str:
        """Prepend and append text to the generated dataset card."""
        content = self.generate(data)
        return prepend_text + '\n' + content + '\n' + postpend_text
            
class DatasetPusher:
    def __init__(self, folder: str, name: str, token: Optional[str] = None, organization: Optional[str] = None, upload_dir: Optional[str] = None):
        """
        Pushes the dataset to the Hugging Face Hub.

        :param folder: Local folder containing dataset files.
        :param name: Name of the dataset repository on the Hugging Face Hub.
        :param token: Hugging Face authentication token (optional, defaults to HF_TOKEN env variable).
        :param organization: Name of the Hugging Face organization (optional).
        :param upload_dir: Directory within the repository where files will be uploaded (optional).
        """
        self.folder = folder
        self.name = name
        self.organization = organization
        self.upload_dir = upload_dir
        self.token = token or os.environ.get("HF_TOKEN")
        if not self.token:
            raise ValueError("Hugging Face token is required. Set it as an environment variable HF_TOKEN or pass it to the constructor.")
        self.api = HfApi()

    def push_to_hub(self, dataset_card: DatasetCard, card_data: Dict[str, Any]):
        current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
        repo_name = f"{self.name}-{current_time}"

        # Determine if it's a personal repo or an organization repo
        if self.organization:
            repo_id = f"{self.organization}/{repo_name}"
        else:
            repo_id = f"{self.api.whoami()['name']}/{repo_name}"

        try:
            self.api.repo_info(repo_id=repo_id, repo_type="dataset")
            print(f"Repository {repo_id} already exists. Pushing to the existing repository.")
        except Exception:
            self.api.create_repo(repo_id=repo_id, repo_type="dataset", private=False)
            print(f"Created new repository: {repo_id}")

        # Setup local repository clone path
        local_dir = f"./temp_{repo_name}"
        repo = Repository(local_dir=local_dir, clone_from=repo_id, repo_type="dataset", use_auth_token=self.token)

        # Walk through the folder and copy files to the local repository
        for root, _, files in os.walk(self.folder):
            for file in files:
                src_path = os.path.join(root, file)

                # Determine the destination path based on the optional upload directory
                if self.upload_dir:
                    dst_path = os.path.join(local_dir, self.upload_dir, os.path.relpath(src_path, self.folder))
                else:
                    dst_path = os.path.join(local_dir, os.path.relpath(src_path, self.folder))

                os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                with open(src_path, "rb") as src, open(dst_path, "wb") as dst:
                    dst.write(src.read())

        # Generate and save the dataset card
        readme_path = os.path.join(local_dir, "README.md")
        dataset_card.save(readme_path, card_data)

        # Commit and push to the Hugging Face Hub
        repo.git_add()
        repo.git_commit("Initial commit")
        repo.git_push()

        print(f"Successfully pushed dataset to https://huggingface.co/datasets/{repo_id}")

        # Clean up the local repository
        repo.delete_repo()
        print("Cleaned up local repository.")

