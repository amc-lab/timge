# Dockerfile
FROM node:20

WORKDIR /app

# Copy package.json and install dependencies
COPY package*.json ./
RUN npm install --production

# Copy the application code
COPY . .

# Build the Next.js application
RUN npm run build

# Expose the default Next.js port
EXPOSE 3000

# Start the application
CMD ["npm", "start"]
