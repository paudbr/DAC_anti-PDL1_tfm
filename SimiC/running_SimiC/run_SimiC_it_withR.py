import argparse
import subprocess as sp
import sys

# Definir la función para manejar argumentos de línea de comandos
def parse_arguments():
    parser = argparse.ArgumentParser(description="Ejecutar script con un directorio específico")
    parser.add_argument("directory", choices=["K25L", "KUV"], help="Directorio a procesar (K25L o KUV)")
    parser.add_argument("subdirectory", choices=["Tumor", "Immune"], help="Subdirectorio a procesar")
    parser.add_argument("treatment", choices=["Combination", "DAC", "PD-L1"], help="Tipo de tratamiento")
    parser.add_argument("lambda1", type=float, help="L1 value")
    parser.add_argument("lambda2",type=float, help="L2 value")
    return parser.parse_args()

args = parse_arguments()
directory = args.directory
subdirectory = args.subdirectory
treatment = args.treatment
l1 = args.lambda1
l2 = args.lambda2

# Asignar el prefijo basado en el directorio
prefix = "KAP25L" if directory == "K25L" else "KAP25L-UV"

def call_simiC(directory, prefix, subdirectory, treatment, l1, l2):
    from simiclasso.clus_regression import simicLASSO_op

    p2df = f'/root/workdir/HDAC_tfm/{directory}/{subdirectory}/inputFiles/{prefix}_controlvs{prefix}_{treatment}_100_1000_subset_matrix.pickle'
    p2assignment = f'/root/workdir/HDAC_tfm/{directory}/{subdirectory}/inputFiles/{prefix}_controlvs{prefix}_{treatment}_phenotype_annotation.txt'
    p2tf = f'/root/workdir/HDAC_tfm/{directory}/{subdirectory}/inputFiles/{prefix}_controlvs{prefix}_{treatment}_100TF_list.pickle'

    p2saved_file = f'/root/workdir/HDAC_tfm/{directory}/{subdirectory}/outputSimic/matrices/Output_{subdirectory}_{directory}_ctrl_{treatment}_L1_{l1}_L2_{l2}_Ws.pickle'
    
    simicLASSO_op(p2df, p2assignment, True, p2tf, p2saved_file, None, 100, 1000, 
                   max_rcd_iter=100000, df_with_label=False,
                   lambda1=l1, lambda2=l2)

    return p2saved_file

# Llamar a la función de SimicLASSO y obtener el archivo generado
output_pickle = call_simiC(directory, prefix, subdirectory, treatment, l1, l2)

print('Filtering weights')
sys.stdout.flush()

# Ejecutar el script de R pasándole el archivo generado
sp.check_call(f'Rscript /root/workdir/HDAC_tfm/scripts/SimiC_bien/Filter_weights_v2.r {output_pickle}', shell=True)

# Generar el nombre del archivo filtrado basado en la nomenclatura usada en el script de R
filtered_pickle = output_pickle.replace('.pickle', '_filtered_BIC.pickle')

def AUC_filter(directory, prefix, subdirectory, treatment, l1, l2, filtered_pickle):
    from simiclasso.weighted_AUC_mat import main_fn

    p2df = f'/root/workdir/HDAC_tfm/{directory}/{subdirectory}/inputFiles/{prefix}_controlvs{prefix}_{treatment}_100_1000_subset_matrix.pickle'
    p2AUC = f'/root/workdir/HDAC_tfm/{directory}/{subdirectory}/outputSimic/matrices/Output_{subdirectory}_{directory}_ctrl_{treatment}_L1_{l1}_L2_{l2}_filtered_AUCs.pickle'

    main_fn(p2df, filtered_pickle, p2AUC)

# Ejecutar la función con la matriz filtrada
AUC_filter(directory, prefix, subdirectory, treatment, l1, l2, filtered_pickle)
