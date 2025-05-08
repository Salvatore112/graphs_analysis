export JAVA_HOME=jav_home

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$JAVA_HOME/jre/lib/amd64/server

export HADOOP_HOME=hadoo_home

export HADOOP_MAPRED_HOME=$HADOOP_HOME

export HADOOP_YARN_HOME=$HADOOP_HOME

export PATH=$PATH:$HADOOP_HOME/bin

export PATH=$PATH:$HADOOP_HOME/sbin

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HADOOP_HOME/lib/native

export CLASSPATH=$CLASSPATH:hadoo_home/hadoop-2.6.1/etc/hadoop

for i in $HADOOP_HOME/share/hadoop/common/lib/*.jar

do

    CLASSPATH=$CLASSPATH:$i

done

for i in $HADOOP_HOME/share/hadoop/hdfs/lib/*.jar

do

    CLASSPATH=$CLASSPATH:$i

done

for i in $HADOOP_HOME/share/hadoop/common/*.jar

do

    CLASSPATH=$CLASSPATH:$i

done

for i in $HADOOP_HOME/share/hadoop/hdfs/*.jar

do

    CLASSPATH=$CLASSPATH:$i

done

export CLASSPATH