import os
from typing import List


class Settings:
    PROJECT_NAME: str = "CFTR-StructMap"
    VERSION: str = "0.1.0"
    DESCRIPTION: str = (
        "Scientific bioinformatics application for reproducible structural-context "
        "analysis of CFTR variants."
    )
    API_V1_STR: str = "/api/v1"

    @property
    def cors_origins(self) -> List[str]:
        raw = os.getenv("CORS_ORIGINS", "http://localhost:3000,http://127.0.0.1:3000")
        return [origin.strip() for origin in raw.split(",") if origin.strip()]


settings = Settings()
