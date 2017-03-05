from msnoise.api import *

from .default import default
from .tomo_table_def import TomoConfig

def main():
    engine = get_engine()
    Session = sessionmaker(bind=engine)
    session = Session()

    try:
        TomoConfig.__table__.create(bind=engine, checkfirst=True)
    except:
        print("Table seems to already exist")

    for name in default.keys():
        try:
            session.add(TomoConfig(name=name, value=default[name][-1]))
            session.commit()
            print("Added %s" % name)
        except:
            print("parameter %s already in database" % name)
            session.rollback()

