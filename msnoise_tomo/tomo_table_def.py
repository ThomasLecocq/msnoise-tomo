# installer for TOMO
from sqlalchemy import Column, String
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class TomoConfig(Base):
    """
    Config Object

    :type name: str
    :param name: The name of the config bit to set.

    :type value: str
    :param value: The value of parameter `name`
    """
    __tablename__ = "tomo-config"
    name = Column(String(255), primary_key=True)
    value = Column(String(255))

    def __init__(self, name, value):
        """"""
        self.name = name
        self.value = value


########################################################################
