from sqlalchemy import create_engine
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import sessionmaker
import logging
import random

logging.basicConfig(level=logging.INFO)

engine = create_engine("sqlite:///flashcards.db", echo=False)
Base = automap_base()
Base.prepare(engine, reflect=True)

Topics = Base.classes.topics
Questions = Base.classes.questions

Session = sessionmaker(bind=engine)
session = Session()

def get_topics():
    try:
        topics = session.query(Topics).all()
        return [topic.name for topic in topics]
    except Exception as e:
        logging.error(f"Error fetching topics: {e}")
        return []
    
def get_topic_question(topic):
    questions = session.query(Questions).filter(Questions.topic_name == topic).all()
    return random.choice([question.question for question in questions])
