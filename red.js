function randomGaussian(mean, stddev) {
    let u = 0, v = 0;
    while (u === 0) u = Math.random(); // Convierte [0,1) a (0,1)
    while (v === 0) v = Math.random();
    const z = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
    return mean + stddev * z;
}



function parent_selection(objects) {
    // Calcula la suma total de "puntos" en todos los objetos
    const totalBajas = objects.reduce((sum, obj) => sum + obj.puntaje, 0);
  
    // Calcula la probabilidad de cada objeto en base a su valor de "puntos"
    const probabilities = objects.map(obj => obj.puntaje / totalBajas);
  
    // Genera un número aleatorio entre 0 y 1
    const random = Math.random();
  
    // Elige el primer objeto con probabilidad ponderada
    let cumulativeProbability = 0;
    let chosenIndex1;
    for (let i = 0; i < objects.length; i++) {
      cumulativeProbability += probabilities[i];
      if (random <= cumulativeProbability) {
        chosenIndex1 = i;
        break;
      }
    }
  
    // Elimina el primer objeto elegido de la lista de objetos
    const remainingObjects = objects.filter((_, index) => index !== chosenIndex1);
  
    // Calcula la suma total de "bajas" en los objetos restantes
    const remainingTotalBajas = remainingObjects.reduce((sum, obj) => sum + obj.puntaje, 0);
  
    // Calcula la probabilidad de cada objeto restante en base a su valor de "bajas"
    const remainingProbabilities = remainingObjects.map(obj => obj.puntaje / remainingTotalBajas);
  
    // Genera un nuevo número aleatorio entre 0 y 1 para el segundo objeto
    const random2 = Math.random();
  
    // Elige el segundo objeto con probabilidad ponderada
    let cumulativeProbability2 = 0;
    let chosenIndex2;
    for (let i = 0; i < remainingObjects.length; i++) {
      cumulativeProbability2 += remainingProbabilities[i];
      if (random2 <= cumulativeProbability2) {
        chosenIndex2 = i;
        break;
      }
    }
  
    // Obtiene los objetos seleccionados
    const chosenObject1 = objects[chosenIndex1];
    const chosenObject2 = remainingObjects[chosenIndex2];
  
    // Retorna los dos objetos seleccionados
    return [chosenObject1, chosenObject2];
  }
  


function printNeuralNetwork(network) {
    console.log('Pesos:');
    for (let i = 0; i < network.weights.length; i++) {
        console.log(`Capa ${i + 1} a Capa ${i + 2}:`);
        for (let j = 0; j < network.weights[i].length; j++) {
            console.log(`Neurona ${j + 1}: ${network.weights[i][j]}`);
        }
    }

    console.log('Sesgos:');
    for (let i = 0; i < network.biases.length; i++) {
        console.log(`Capa ${i + 2}: ${network.biases[i]}`);
    }
}

function printGenome(genome) {
    console.log('Genoma:');
    console.log('Layers:', genome.layers);
    console.log('Número de genes:', genome.num_genes);
    console.log('Genes:');
    for (let i = 0; i < genome.genes.length; i++) {
        const gene = genome.genes[i];
        console.log(`Gen ${i + 1}:`);
        console.log('  Neurona ID1:', gene.neuron_id1);
        console.log('  Neurona ID2:', gene.neuron_id2);
        console.log('  Capa ID1:', gene.layer_id1);
        console.log('  Capa ID2:', gene.layer_id2);
        console.log('  Habilitado:', gene.enabled);
        console.log('  Peso:', gene.weight);
        console.log('  Sesgo 1:', gene.bias1);
        console.log('  Sesgo 2:', gene.bias2);
    }
}
function create_neural_network(genome_r, list_networks) {
    var genome = Object.create(Object.getPrototypeOf(genome_r));
    Object.assign(genome, genome_r);
    var layers = genome.layers;
    var genes = genome.genes;
    var network;
    var weights, biases;
  
    if (list_networks === undefined) {
      network = new NeuralNetwork(layers, genome, true);
      weights = network.get_weights();
      biases = network.get_biases();
    } else {
      network = list_networks.shift();
      network.genome = genome;
      weights = [];
      biases = [];
  
      for (var i = 0; i < network.layer_sizes.length - 1; i++) {
        weights.push(new Array(network.layer_sizes[i + 1]).fill(0).map(function () {
          return new Array(network.layer_sizes[i]).fill(0);
        }));
        biases.push(new Array(network.layer_sizes[i + 1]).fill(0));
      }
    }
  
    for (var _i = 0; _i < genes.length; _i++) {
      var gene = genes[_i];
  
      if (gene.enabled) {
        var layer_idx = gene.layer_id1;
        var neuron_idx1 = gene.neuron_id1;
        var neuron_idx2 = gene.neuron_id2;
        var weight = gene.weight;
        var bias1 = gene.bias1;
        var bias2 = gene.bias2;
  
        if (layer_idx < layers.length - 1) {
          weights[layer_idx][neuron_idx2][neuron_idx1] = weight;
  
          if (layer_idx !== 0) {
            biases[layer_idx - 1][neuron_idx1] = bias1;
          }
  
          biases[layer_idx][neuron_idx2] = bias2;
        }
      }
    }
  
    network.set_weights(weights);
    network.set_biases(biases);
    return network;
}

// Clase Gen
class Gene {
    constructor(neuron_id1, neuron_id2, layer_id1, layer_id2, bias1 = 0, bias2 = 0, weight = 0, random = false) {
        this.neuron_id1 = neuron_id1;
        this.neuron_id2 = neuron_id2;
        this.layer_id1 = layer_id1;
        this.layer_id2 = layer_id2;

        this.enabled = true;

        if (random) {
            this.weight = randomGaussian(0,1)
            this.bias1 = randomGaussian(0,1)
            this.bias2 = randomGaussian(0,1)
        } else {
            this.weight = weight;
            this.bias1 = bias1;
            this.bias2 = bias2;
        }
    }
}


// Clase Genoma
class Genome {
    constructor(layers, num_genes, rand, debug=false) {
        this.layers = layers;
        this.num_genes = num_genes;
        this.genes = this.generate_genes(rand);

        // Probabilidad de cambiar los parámetros como los sesgos y el peso de la conexión
        this.mutation_rate_parameter = 0;
        // Probabilidad de cambiar uno de los nodos al que se conecta
        this.mutation_rate_connection = 0;
        // Probabilidad de deshabilitar una conexión habilitada
        this.mutation_rate_disabled = 0;
        // Probabilidad de habilitar una conexión deshabilitada
        this.mutation_rate_enabled = 0;
        //probabilidad de que se añada un gen
        this.mutation_rate_add = 0
        // Probabilidad de que se elimine un gen
        this.mutation_rate_remove = 0

        // Probabilidad de cruce genetico
        this.cruce = 0

        this.debug = debug
    }

    generate_genes(rand) {
        const genes = [];
        const num_layers = this.layers.length;

        for (let i = 0; i < this.num_genes; i++) {
            const capaI = Math.floor(Math.random() * (num_layers - 1));
            const capaF = capaI + 1;
            let neuronasAnt = 0;
            for (let j = 0; j < capaI; j++) {
                neuronasAnt += this.layers[j];
            }

            const gene = new Gene(
                Math.floor(Math.random() * this.layers[capaI]),
                Math.floor(Math.random() * this.layers[capaF]),
                capaI,
                capaF,
                0,
                0,
                0,
                rand
            );
            genes.push(gene);
        }

        return genes;
    }

    crossover(other_genome) {
        const new_genome = new Genome(this.layers, {... this.num_genes}, false, this.debug);
        const new_genome2 = new Genome(this.layers, {... this.num_genes}, false,this.debug);

        new_genome.mutation_rate_parameter = this.mutation_rate_parameter
        new_genome.mutation_rate_connection = this.mutation_rate_connection
        new_genome.mutation_rate_disabled = this.mutation_rate_disabled
        new_genome.mutation_rate_enabled = this.mutation_rate_enabled
        new_genome.mutation_rate_remove = this.mutation_rate_remove 
        new_genome.mutation_rate_add = this.mutation_rate_add 
        
        new_genome2.mutation_rate_parameter = this.mutation_rate_parameter
        new_genome2.mutation_rate_connection = this.mutation_rate_connection
        new_genome2.mutation_rate_disabled = this.mutation_rate_disabled
        new_genome2.mutation_rate_enabled = this.mutation_rate_enabled
        new_genome2.mutation_rate_remove = this.mutation_rate_remove 
        new_genome2.mutation_rate_add = this.mutation_rate_add 

        if (Math.random() < this.cruce){
            const crossover_point = Math.floor(Math.random() * Math.min(this.genes.length, other_genome.genes.length));

            new_genome.genes = other_genome.genes.slice(0, crossover_point).concat(this.genes.slice(crossover_point));
            new_genome2.genes = this.genes.slice(0, crossover_point).concat(other_genome.genes.slice(crossover_point));
        } else {
            new_genome.genes = this.genes.slice();
            new_genome2.genes = other_genome.genes.slice();
        }


        // console.log(new_genome.genes.length)


        

        return [new_genome, new_genome2];
    }
    crossoverOnePoint(other_genome) {
        const new_genome = new Genome(this.layers, this.num_genes, false, this.debug);
        const new_genome2 = new Genome(this.layers, this.num_genes, false, this.debug);
    
        const crossover_point = Math.floor(Math.random() * Math.min(this.genes.length, other_genome.genes.length));
    
        new_genome.genes = this.genes.slice(0, crossover_point).concat(other_genome.genes.slice(crossover_point));
        new_genome2.genes = other_genome.genes.slice(0, crossover_point).concat(this.genes.slice(crossover_point));
    
        return [new_genome, new_genome2];
    }
    crossoverTwoPoints(other_genome) {
        const new_genome = new Genome(this.layers, this.num_genes, false, this.debug);
        const new_genome2 = new Genome(this.layers, this.num_genes, false, this.debug);
    
        const crossover_point1 = Math.floor(Math.random() * Math.min(this.genes.length, other_genome.genes.length));
        const crossover_point2 = Math.floor(Math.random() * Math.min(this.genes.length, other_genome.genes.length));
    
        const start_point = Math.min(crossover_point1, crossover_point2);
        const end_point = Math.max(crossover_point1, crossover_point2);
    
        new_genome.genes = this.genes.slice(0, start_point).concat(other_genome.genes.slice(start_point, end_point)).concat(this.genes.slice(end_point));
        new_genome2.genes = other_genome.genes.slice(0, start_point).concat(this.genes.slice(start_point, end_point)).concat(other_genome.genes.slice(end_point));
    
        return [new_genome, new_genome2];
    }
    crossoverUniform(other_genome) {
        const new_genome = new Genome(this.layers, this.num_genes, false, this.debug);
        const new_genome2 = new Genome(this.layers, this.num_genes, false, this.debug);
    
        for (let i = 0; i < this.genes.length; i++) {
            if (Math.random() < 0.5) {
                new_genome.genes.push(this.genes[i]);
                new_genome2.genes.push(other_genome.genes[i]);
            } else {
                new_genome.genes.push(other_genome.genes[i]);
                new_genome2.genes.push(this.genes[i]);
            }
        }
    
        return [new_genome, new_genome2];
    }
    mutate() {
        let p = 0; // Mutaciones de parámetros (pesos, sesgos)
        let c = 0; // Mutaciones de conexiones (cambiar los nodos conectados)
        let dis = 0; // Mutaciones de deshabilitar una conexión
        let en = 0; // Mutaciones de habilitar una conexión
        let add = 0; // Mutaciones de agregar un nuevo gen
        let rem = 0; // Mutaciones de eliminar un gen
    
        for (const gene of this.genes) {
            if (Math.random() < this.mutation_rate_parameter) {
                gene.weight += randomGaussian(0, 1);
                gene.bias1 += randomGaussian(0, 1);
                gene.bias2 += randomGaussian(0, 1);
                p++;
            }
    
            if (Math.random() < this.mutation_rate_connection) {
                c++;
    
                if (Math.random() < 0.5) {
                    gene.neuron_id2 = Math.floor(Math.random() * this.layers[gene.layer_id1 + 1]);4
                } else {
                    gene.neuron_id1 = Math.floor(Math.random() * this.layers[gene.layer_id1]);
                }
            }
        }
    
        if (Math.random() < this.mutation_rate_disabled) {
            dis++;
    
            for (const g of this.genes) {
                if (g.enabled) {
                    g.enabled = false;
                    break;
                }
            }
        }
    
        if (Math.random() < this.mutation_rate_enabled) {
            en++;
    
            for (const g of this.genes) {
                if (!g.enabled) {
                    g.enabled = true;
                    break;
                }
            }
        }
    
        if (Math.random() < this.mutation_rate_add) {
            add++;
    
            const capaI = Math.floor(Math.random() * (this.layers.length - 1));
            const capaF = capaI + 1;
            const gene = new Gene(
                Math.floor(Math.random() * this.layers[capaI]),
                Math.floor(Math.random() * this.layers[capaF]),
                capaI,
                capaF,
                randomGaussian(0, 1),
                randomGaussian(0, 1),
                randomGaussian(0, 1),
                true
            );
    
            this.genes.push(gene);
            this.num_genes+=1
        }
    
        if (Math.random() < this.mutation_rate_remove && this.genes.length > 1) {
            rem++;
    
            const index = Math.floor(Math.random() * this.genes.length);
            this.genes.splice(index, 1);
            this.num_genes-=1
        }
    
        if (this.debug) {
            console.log("Mutaciones:\n");
            console.log("Parámetro =", p);
            console.log("Conexión =", c);
            console.log("Deshabilitar =", dis);
            console.log("Habilitar =", en);
            console.log("Agregar =", add);
            console.log("Eliminar =", rem);
            console.log("");
        }
    }
}


// Red Neuronal
class NeuralNetwork {
    constructor(layer_sizes, genome, zeros = false) {
        this.num_layers = layer_sizes.length;
        this.layer_sizes = layer_sizes;
        if (zeros) {
            this.weights = this.generateWeightArray(layer_sizes);
            this.biases = layer_sizes.slice(1).map(i => Array(i).fill().map(() => 0));
        } else {
            this.weights = this.generateWeightArray(layer_sizes);
            this.biases = layer_sizes.slice(1).map(i => Array(i).fill().map(() => Math.random()));
        }
        this.genome = genome;

        this.capas_ocultas = activation_functions.sigmoid;
        this.capa_salida = activation_functions.sigmoid;
    }
    generateWeightArray(layer_sizes,rand=false, mean = 0, stddev = 1, isBias = false) {
        const numLayers = layer_sizes.length;
        const weightArray = [];
    
        for (let i = 1; i < numLayers; i++) {
            const prevLayerSize = layer_sizes[i - 1];
            const currentLayerSize = layer_sizes[i];
            const layerWeights = [];
    
            for (let j = 0; j < currentLayerSize; j++) {
                const neuronWeights = [];
                let weight
                for (let k = 0; k < prevLayerSize; k++) {
                    if (rand){
                        weight = randomGaussian(mean, stddev);
                    }else{
                        weight = 0;
                    }
                    neuronWeights.push(weight);
                }
    
                layerWeights.push(neuronWeights);
            }
    
            weightArray.push(layerWeights);
        }
    
        if (isBias) {
            return weightArray.map(layer => layer.map(neuronWeights => mean));
        }
    
        return weightArray;
    }
    
    forward(inputs) {
        let activations = inputs;

        for (let i = 0; i < this.num_layers -1; i++) {
            const layer_weights = this.weights[i];
            const layer_biases = this.biases[i];
            const layer_activation = [];
            for (let j = 0; j < layer_weights.length; j++) {
                let weighted_sum = 0;

                for (let k = 0; k < layer_weights[j].length; k++) {
                    weighted_sum += layer_weights[j][k] * activations[k];
                }
                weighted_sum += layer_biases[j];
                if (i<this.num_layers-2){
                    layer_activation.push(this.capas_ocultas(weighted_sum));
                }else{
                    let r = this.capa_salida(weighted_sum)
                    layer_activation.push(r);
                }
            }
            // print("capa")

            activations = {...layer_activation};
        }
        // print()
        return activations;
    }

    set_weights(weights) {
        if (weights.length !== this.weights.length) {
            throw new Error("Number of weight matrices does not match.");
        }
        this.weights = weights;
    }

    set_weights_wI(layer_index, weights) {
        if (weights.length !== this.weights[layer_index].length) {
            throw new Error("Number of weight matrices does not match.");
        }
        this.weights[layer_index] = weights;
    }

    set_biases(biases) {
        if (biases.length !== this.biases.length) {
            throw new Error("Number of bias vectors does not match.");
        }
        this.biases = biases;
    }

    get_weights() {
        return this.weights;
    }

    get_biases() {
        return this.biases;
    }

    apply_mutations(mutation_rate) {
        const num_weights = this.weights.length;
        const num_biases = this.biases.length;

        for (let i = 0; i < num_weights; i++) {
            if (Math.random() < mutation_rate) {
                this.weights[i] = matrixAdd(this.weights[i], matrixRandomUniform(this.weights[i].length, this.weights[i][0].length, -1, 1));
            }
        }

        for (let i = 0; i < num_biases; i++) {
            if (Math.random() < mutation_rate) {
                this.biases[i] = matrixAdd(this.biases[i], matrixRandomUniform(this.biases[i].length, this.biases[i][0].length, -1, 1));
            }
        }
    }
}

// Helper function: Dot product of two matrices
function matrixDot(a, b) {
    const result = Array(a.length).fill().map(() => Array(b[0].length).fill(0));
    for (let i = 0; i < a.length; i++) {
        for (let j = 0; j < b[0].length; j++) {
            for (let k = 0; k < a[0].length; k++) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

// Helper function: Add two matrices element-wise
function matrixAdd(a, b) {
    return a.map((row, i) => row.map((cell, j) => cell + b[i][j]));
}

// Helper function: Create a matrix with random values from a uniform distribution
function matrixRandomUniform(rows, cols, min, max) {
    return Array(rows).fill().map(() => Array(cols).fill().map(() => Math.random() * (max - min) + min));
}

// Activation Functions
const activation_functions = {
    sigmoid: x => 1 / (1 + Math.exp(-x)),
    relu: x => Math.max(0, x),
    tanh: x => Math.tanh(x)
};


class Rect {
    constructor(pos, width, height) {
        this.pos = pos
        this.width = width;
        this.height = height;
        
        this.activo = true
    }
    
    // Método para comprobar si este rectángulo colisiona con otro rectángulo
    isCollidingWith(otherRect) {
        // Calcula las coordenadas de los bordes del rectángulo actual
        const left = this.pos[0] - this.width / 2;
        const right = this.pos[0] + this.width / 2;
        const top = this.pos[1] - this.height / 2;
        const bottom = this.pos[1] + this.height / 2;
        
        // Calcula las coordenadas de los bordes del otro rectángulo
        const otherLeft = otherRect.pos[0] - otherRect.width / 2;
        const otherRight = otherRect.pos[0] + otherRect.width / 2;
        const otherTop = otherRect.pos[1] - otherRect.height / 2;
        const otherBottom = otherRect.pos[1] + otherRect.height / 2;
        
        // Comprueba si hay superposición entre los rectángulos en los ejes x e y
        if (left < otherRight && right > otherLeft && top < otherBottom && bottom > otherTop && this.activo && otherRect.activo) {
            return true; // Hay colisión
        } else {
            return false; // No hay colisión
        }
    }
}

function getRandomNumber(min, max) {
    const random = Math.random();
    const randomNumber = Math.floor(random * (max - min + 1)) + min;
    
    return randomNumber;
}