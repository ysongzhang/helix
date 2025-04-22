NPC=$1

IPs=(
172.16.0.2
172.16.0.3
172.16.0.4
172.16.0.5
172.16.0.6
172.16.0.7
172.16.0.8
172.16.0.9
172.16.0.10
172.16.0.11
172.16.0.12
)

file_name=Inference/IP_HOSTS/IP_$NPC
if [ $NPC = 11 ]; then
    start_port=5000
    for((i=0;i<11;i++)); do
        echo ${IPs[i]}:$start_port >> $file_name
        start_port=$(($start_port + 1))
        
    done

elif [ $NPC = 21 ]; then
    start_port=5000
    for((i=0;i<10;i++)); do
        echo ${IPs[i]}:$start_port >> $file_name
        start_port=$(($start_port + 1))
        echo ${IPs[i]}:$start_port >> $file_name
        start_port=$(($start_port + 1))
        
    done
    echo ${IPs[10]}:$start_port >> $file_name

elif [ $NPC = 31 ]; then
    start_port=5000
    for((i=0;i<10;i++)); do
        echo ${IPs[i]}:$start_port >> $file_name
        start_port=$(($start_port + 1))
        echo ${IPs[i]}:$start_port >> $file_name
        start_port=$(($start_port + 1))
        echo ${IPs[i]}:$start_port >> $file_name
        start_port=$(($start_port + 1))
        
    done
    echo ${IPs[10]}:$start_port >> $file_name

elif [ $NPC = 63 ]; then
    start_port=5000
    for((i=0;i<8;i++)); do
        for((j=0;j<6;j++)); do
            echo ${IPs[i]}:$start_port >> $file_name
            start_port=$(($start_port + 1))
        done
        
    done

    for((i=8;i<11;i++)); do
        for((j=0;j<5;j++)); do
            echo ${IPs[i]}:$start_port >> $file_name
            start_port=$(($start_port + 1))
        done
        
    done
    
elif [ $NPC = 127 ]; then
    start_port=5000
    for((i=0;i<6;i++)); do
        for((j=0;j<12;j++)); do
            echo ${IPs[i]}:$start_port >> $file_name
            start_port=$(($start_port + 1))
        done
        
    done

    for((i=6;i<11;i++)); do
        for((j=0;j<11;j++)); do
            echo ${IPs[i]}:$start_port >> $file_name
            start_port=$(($start_port + 1))
        done
        
    done

elif [ $NPC = 255 ]; then
    start_port=5000
    for((i=0;i<2;i++)); do
        for((j=0;j<24;j++)); do
            echo ${IPs[i]}:$start_port >> $file_name
            start_port=$(($start_port + 1))
        done
        
    done

    for((i=2;i<11;i++)); do
        for((j=0;j<23;j++)); do
            echo ${IPs[i]}:$start_port >> $file_name
            start_port=$(($start_port + 1))
        done
        
    done
    
fi