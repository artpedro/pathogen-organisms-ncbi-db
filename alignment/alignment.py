# pseudo-code

# verificar se o arquivo de entrada é válido

# definir a classe ranking
    
    # atributos: ranking = [(patógeno,assembly),...];

    # métodos: add_entry(asm_acc,patógeno,porcentagem): adiciona o resultado do alinhamento no ranking;
    #          add_batch([(patógeno, assembly, porcentagem)]): adiciona o resultado do alinhamento de um patógeno inteiro;
    #          retrieve_entry(position): devolve qual a tupla nessa posição 0 - 99;

# função alinhar(index,input)
    # recebe o index e o input e realiza o alinhamento
    # salva o output no temp.log
    # devolve a porcentagem do alinhamento


# para cada patógeno
    # salvar o caminho pro index
    # listar todos o index presentes
    
    # declarar o batch
    
    {
    
        # para cada index
        # armazenar o asm
        # realizar o alinhamento e armazenar o stdout e stderr em temp.log
        # excluir .sam
        # adicionar o registro em batch
    
    } descobrir como paralelizar essa parte

    # adicionar o batch ao ranking

# verificar qual são os 10 melhores alinhamentos
# refazer os 10 alinhamentos e salvar o .sam

