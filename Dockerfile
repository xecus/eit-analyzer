FROM gcc:4.9
COPY . /app
WORKDIR /app

# Write in progress...

CMD ["gcc -v"]
