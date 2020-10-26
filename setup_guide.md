# Reactome Setup Guide

The following are the instructions to get Reactome running on Ubuntu 18.04. For other OS, please refer to the official Neo4j documentation.

### 1. Install Java 8

Neo4j specifically requires Java 8 to be installed. If you don't have Java installed, you can install it with the command below:
```bash
$ sudo apt install openjdk-8-jdk
```
Otherwise validate that you have Java 8 installed:
```bash
$ java -version
openjdk version "1.8.0_181"
```

### 2. Install Neo4j Community Edition

[Linux Installation - Neo4j Reference](https://neo4j.com/docs/operations-manual/current/installation/linux/debian/?_ga=2.249168388.2041192375.1507250087-893468657.1507250087).

Run the following commands to install Neo4j Community Edition 3.4.6:
```bash
$ wget -O - https://debian.neo4j.org/neotechnology.gpg.key | sudo apt-key add -
$ echo 'deb https://debian.neo4j.org/repo stable/' | sudo tee -a /etc/apt/sources.list.d/neo4j.list
$ sudo apt-get update
$ sudo apt-get install neo4j=1:3.4.6
```
Later version of Neo4j can also be used, as long as it is version 3 (version 4 seems to have problems with Reactome database).

Verify that Neo4j is running:
```bash
$ sudo service neo4j status
```
From the status above, you can see that $NEO4J_HOME is located at `/var/lib/neo4j`. 
If Neo4j is not running, start it:
```bash
$ sudo service neo4j start
```

Once it's running, [set the initial password](https://stackoverflow.com/questions/47530154/neo4j-command-failed-initial-password-was-not-set-because-live-neo4j-users-wer) to whatever you prefer.
```bash
$ curl -H "Content-Type: application/json" -X POST -d '{"password":"WHATEVER THE PASSWORD IS"}' -u neo4j:neo4j http://localhost:7474/user/neo4j/password
```

### 3. Install Reactome database

See https://reactome.org/dev/graph-database

Download the Reactome database. Extract and move it to `$NEO4J_HOME/data/databases`.
```bash
$ wget https://reactome.org/download/current/reactome.graphdb.tgz
$ tar xvzf reactome.graphdb.tgz
$ sudo mv graph.db /var/lib/neo4j/data/databases
$ chown -R neo4j:neo4j /var/lib/neo4j/data/databases/graph.db
```
Edit the config file at either `$NEO4J_HOME/conf/neo4j.conf` or `/etc/neo4j/neo4j.conf`. 
Change ```dbms.active_database``` to ```dbms.active_database=graph.db``` if necessary.

Check that the neo4j service is running with the following command. If it is not running, start it.
```bash
$ sudo service neo4j status
```

For graph database connection in PALS, be sure to set the following environmental variables:
- `NEO4J_SERVER`: your Neo4j server (default: bolt://localhost:7687)
- `NEO4J_USER`: your Neo4j user name (default: neo4j)
- `NEO4J_PASSWORD`: your Neo4j password (default: neo4j)