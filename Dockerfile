FROM python:3.8-slim-buster

WORKDIR /app

RUN python -m pip install pipenv

COPY ./ ./

RUN pipenv install

EXPOSE 3003

CMD ["pipenv", "run", "python", "-u", "backend/search.py"]