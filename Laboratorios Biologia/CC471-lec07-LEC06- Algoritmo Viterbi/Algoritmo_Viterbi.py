def AlgoritmoViterbi(secuencia):
    """Recibe una secuencia de salida y retona la secuencia de estados del HMM con la 
       probabilidad más alta, dada la secuencia de simbolos de salida."""
    
    #=====================================================================================
    # Las regiones son I(Island) y G(Genome)
    region = ['I', 'G']   
    nucleotidos = ['C', 'G', 'A', 'T']
    #Creamos la probabiliada de emision, transicion, iniciacion
    p_e = [[0.4, 0.4, 0.1, 0.1], 
           [0.1, 0.1, 0.4, 0.4]]
    p_t = [[0.8, 0.2], 
           [0.1, 0.9]]
    p_i = [0.1, 0.9]   

    # Creamos un diccionario por cada matriz
    prob_e = {}  ;  prob_t = {}   ;   prob_i = {}
    i = 0
    for reg in region:
        # Diccionario para la probabilidad de emision y transmicion
        dic_e = {}   ;  dic_t = {} 
        j = 0
        for nuc in nucleotidos:
            dic_e.update({nuc : p_e[i][j]}) ; j += 1
        j = 0
        for new_reg in region:
            dic_t.update({new_reg : p_t[i][j]}) ; j += 1
        prob_e.update({reg : dic_e})
        prob_t.update({reg : dic_t})
        prob_i.update({reg : p_i[i]})
        i += 1
    
    #=====================================================================================
    # Creamos una lista para cada estado, esta lista contiene el nucleotido de salida y
    # un diccionario que contiene diccionarios por cada region G e I, el dict de cada region 
    # contiene la probabilidad hasta ese instante(Programacion dinamica) y el origen de este estado(Backtraking)  
    # S1 = ['nucleotido', dict{}] , dict_reg = {'region': dict_to_reg }  ,  dict_to_reg = {probabilidad, origen}
    #=====================================================================================
    
    estado = [] # Diccionario de todos los estados 

    for i in range(len(secuencia)):
        s = {}   ;  dict_reg = {}; 
        # Para el estado inicial 
        if i == 0:
            # Completamos las probabilidades de todas las regiones del estado inicial
            for reg in region:
                dict_to_reg = {}     # Diccionario por cada region
                nucleotido_inicial = secuencia[i]
                p = prob_i[reg] * prob_e[reg][nucleotido_inicial]
                # Agregamos al diccionario de cada region
                dict_to_reg.update({'Prob' : p})
                dict_to_reg.update({'Origen' : ' '})
                # Agregamos al diccionario general va contener todas las regiones
                dict_reg.update({reg : dict_to_reg})
            # Agregamos el nucleotido y el diccionario total a la lista del estado
            s.update({'Salida' : secuencia[i]});   
            s.update({'Datos' : dict_reg})
            # Agregamos es estado s a la lista de todos los estados
            #estado.update({'S_'+str(i): s})
            estado.append(s)   

        # Para los demas estados
        else:
            # Hallamos la probabilidad para cada region{I, G} del nuevo estado
            for new_region in region:
                dict_to_reg = {}     # Diccionario por cada region
                p_final = 0 ; origen_final = ' ' 
                # Necesitamos la maxima prob de todas las regiones de estado anterior 
                for old_region in region:
                    p_s_anterior = estado[i-1]['Datos'][old_region]['Prob']
                    p_transicion = prob_t[old_region][new_region]
                    nucleotido = secuencia[i]
                    p_emision = prob_e[new_region][nucleotido]
                    p = p_s_anterior * p_transicion * p_emision
                    if p > p_final:
                        p_final = p   ;   origen_final = old_region

                # Agregamos al diccionario de cada region
                dict_to_reg.update({'Prob' : p_final})
                dict_to_reg.update({'Origen' : origen_final})
                # Agregamos al diccionario general va contener todas las regiones
                dict_reg.update({new_region : dict_to_reg})

            # Agregamos el nucleotido y el diccionario total a la lista del estado
            s.update({'Salida' : secuencia[i]});   
            s.update({'Datos' : dict_reg})
            # Agregamos es estado s a la lista de todos los estados
            #estado.update({'S_'+str(i): s})
            estado.append(s)
    #Mostrar los estados
    #for est in estado:
    #    print(est)
    
    #=====================================================================================
    
    # Backtraking: encontrar el mejor camino
    # En el ultimo estado se encuentra la mayor probabilidad 

    path = []

    # Buscar el maximo en todas las regiones del ultimo estado
    prob_max = 0      ;   reg_max_prob = ' '
    for reg in region:
        p = estado[-1]['Datos'][reg]['Prob']
        if p > prob_max:
            prob_max = p   ; reg_max_prob = reg

    pivote  = reg_max_prob
    path.append(pivote)
    #Recorremos desde el ultimo a el primero
    for i in reversed(range(len(estado))):
    # Con el valor de path podemos saber donde nos encontramos
    # en el backtraking
        pivote = estado[i]['Datos'][pivote]['Origen']
        path.append(pivote)

    path = path[:-1]
    #Revertir el path
    path_final = []
    for i in reversed(path):
        path_final.append(i)
    return path


if __name__ == '__main__' :
    print("Prueba con secuencia del PPT de la clase") 
    secuencia = ['G', 'T', 'G', 'C', 'C', 'T', 'A']
    path_final = AlgoritmoViterbi(secuencia)
    print("Camino obtimo:",path_final)
    print("Secuencia_1:  ", secuencia,"\n")
    
    secuencia1 = ['T', 'A', 'C', 'G', 'A']
    path_final = AlgoritmoViterbi(secuencia1)
    print("Camino obtimo:",path_final)
    print("Secuencia_2:  ", secuencia1,"\n")
    
    secuencia2 = ['A', 'T', 'C', 'C', 'A', 'T', 'G', 'C', 'G']
    path_final = AlgoritmoViterbi(secuencia2)
    print("Camino obtimo:",path_final)
    print("Secuencia_3:  ", secuencia2)
    