"""Database session and engine management."""

import os
from pathlib import Path
from typing import Generator, Optional

from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from .models import Base

DEFAULT_DB_PATH = Path("data/processed/cftr_structmap.db")


def get_database_url(db_input: Optional[str] = None) -> str:
    """Resolve database URL from parameter, environment, or default path.
    
    Handles raw filesystem paths (converts to sqlite:///), in-memory specifiers,
    or pre-formatted SQLAlchemy connection URLs.
    """
    if db_input:
        if db_input == ":memory:":
            return "sqlite:///:memory:"
        if "://" in db_input:
            return db_input
        # Plain filesystem path provided
        return f"sqlite:///{Path(db_input).resolve()}"

    env_url = os.getenv("DATABASE_URL")
    if env_url:
        return env_url

    DEFAULT_DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    return f"sqlite:///{DEFAULT_DB_PATH.resolve()}"


def get_engine(db_input: Optional[str] = None):
    """Create SQLAlchemy engine with appropriate SQLite connect_args."""
    url = get_database_url(db_input)
    connect_args = {"check_same_thread": False} if url.startswith("sqlite") else {}
    return create_engine(url, connect_args=connect_args)


def init_db(engine=None):
    """Create all database tables."""
    eng = engine or get_engine()
    Base.metadata.create_all(bind=eng)
    return eng


def get_session_factory(engine=None):
    """Get a configured sessionmaker."""
    eng = engine or get_engine()
    return sessionmaker(autocommit=False, autoflush=False, bind=eng)


def get_db(db_input: Optional[str] = None) -> Generator[Session, None, None]:
    """Dependency / context generator for database sessions."""
    engine = get_engine(db_input)
    session_factory = get_session_factory(engine)
    session = session_factory()
    try:
        yield session
    finally:
        session.close()
