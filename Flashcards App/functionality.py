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
    
def get_topic_question(topic, current_question, include_known):
    query = session.query(Questions).filter(Questions.topic_name == topic)
    
    if not include_known:
        query = query.filter(Questions.is_known == 0)

    questions = query.all()

    if not questions:
        return None
    
    filtered_questions = [question.question for question in questions if question.question != current_question]

    return random.choice(filtered_questions) if filtered_questions else random.choice([question.question for question in questions])

def get_question_answer(question_text):
    question = session.query(Questions).filter(Questions.question == question_text).first()

    if question:
        return question.answer
    else:
        return None
    
def set_question_known(is_known, question_text):
    question = session.query(Questions).filter(Questions.question == question_text).first()
    
    if question:
        question.is_known = is_known
        session.commit()

def create_topic(topic_name, warning_label):
    if len(topic_name) < 1:
        warning_label.setText("Name cant be empty!")
        return
    is_existing = session.query(Topics).filter_by(name=topic_name).first()
    if is_existing:
        warning_label.setText("Topic like this already exists!")
        return
    new_topic = Topics(name = topic_name)
    session.add(new_topic)
    session.commit()
    return warning_label.setText("Set created!")
