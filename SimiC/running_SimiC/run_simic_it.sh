#!/bin/bash
#!/bin/bash
# --------------------------------------------------------------------
# Title: run_simic_it.sh
# Author: Paula de Blas Rioja
#(INITIAL VERSION)
# Description: Run SimiC analysis for all specified directory/treatment combinations.
# Usage: ./run_simic_it_all.sh lambda1 lambda2
# --------------------------------------------------------------------
# Instructions:
# 1. Load required Python environment and dependencies.
# 2. Set run paths and parameters (lambda1, lambda2, directories, subdirectories).
# 3. Run the SimiC script for each combination of directory and subdirectory.
# 4. Save the output weight matrix and logs.
# --------------------------------------------------------------------


# Verificar si se proporcionan exactamente 3 argumentos
if [ "$#" -ne 3 ]; then
    echo "Uso: $0 {K25L|KUV} lambda1 lambda2"
    exit 1
fi

directory=$1  # Capturar el primer argumento
lambda1=$2    # Capturar el segundo argumento
lambda2=$3    # Capturar el tercer argumento

# Crear directorio de logs si no existe
mkdir -p ~/workdir/HDAC_tfm/logs

# Definir listas de subdirectorios y tratamientos
subdirectories=("Tumor" "Immune")
treatments=("DAC" "PD-L1" "Combination")

# Iterar sobre todas las combinaciones y ejecutar el script Python
for subdirectory in "${subdirectories[@]}"; do
    for treatment in "${treatments[@]}"; do
        log_file="$HOME/workdir/HDAC_tfm/logs/simic_run_${directory}_${subdirectory}_ctrl_${treatment}_L1_${lambda1}_L2_${lambda2}.log"
        echo "Ejecutando: prefix=$directory, subdirectory=$subdirectory, treatment=$treatment"

        # Ejecutar con nohup y obtener el PID
        nohup python3 ~/workdir/HDAC_tfm/scripts/SimiC_bien/run_SimiC_it_withR.py "$directory" "$subdirectory" "$treatment" "$lambda1" "$lambda2"> "$log_file" 2>&1 &
        pid=$!

        # Esperar a que el proceso termine antes de continuar con la siguiente combinación
        wait $pid



        echo "Finalizado: prefix=$directory, subdirectory=$subdirectory, treatment=$treatment"
    done
done

