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
        return self.template.format(**data)

    def save(self, output_path: str, data: Dict[str, Any]):
        content = self.generate(data)
        with open(output_path, 'w') as f:
            f.write(content)
            
class DatasetPusher:
    def __init__(self, folder: str, name: str, token: Optional[str] = None):
        """
        Made to push our first dataset to the hub!

        Usage example:
        pusher = DatasetPusher(folder="path/to/your/dataset/folder", name="my-dataset")
        pusher.push_to_hub()
        """
        self.folder = folder
        self.name = name
        self.token = token or os.environ.get("HF_TOKEN")
        if not self.token:
            raise ValueError("Hugging Face token is required. Set it as an environment variable HF_TOKEN or pass it to the constructor.")
        self.api = HfApi()

    def push_to_hub(self, dataset_card: DatasetCard, card_data: Dict[str, Any]):
        current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
        repo_name = f"{self.name}-{current_time}"

        try:
            self.api.repo_info(repo_id=repo_name, repo_type="dataset")
            print(f"Repository {repo_name} already exists. Pushing to existing repository.")
        except Exception:
            self.api.create_repo(repo_id=repo_name, repo_type="dataset", private=False)
            print(f"Created new repository: {repo_name}")

        local_dir = f"./temp_{repo_name}"
        repo = Repository(local_dir=local_dir, clone_from=f"{self.api.whoami()['name']}/{repo_name}", repo_type="dataset", use_auth_token=self.token)

        for root, _, files in os.walk(self.folder):
            for file in files:
                src_path = os.path.join(root, file)
                dst_path = os.path.join(local_dir, os.path.relpath(src_path, self.folder))
                os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                with open(src_path, "rb") as src, open(dst_path, "wb") as dst:
                    dst.write(src.read())

        # Generate and save the dataset card
        readme_path = os.path.join(local_dir, "README.md")
        dataset_card.save(readme_path, card_data)

        repo.git_add()
        repo.git_commit("Initial commit")
        repo.git_push()

        print(f"Successfully pushed dataset to https://huggingface.co/datasets/{self.api.whoami()['name']}/{repo_name}")

        repo.delete_repo()
        print("Cleaned up local repository.")

