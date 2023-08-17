let canvas, ctx, birds, birds_Eliminados
let fondo = []
let suelo = []
let rectSuelo
let teclas = {};
let teclasUp = {};
radians = (Math.PI / 180)
let reduccion_giro = 20
let suelo_x = 0
let velocidad_juego = 0
let tubos = []
let distancia_regeneracion = 230
let puntaje = 0
let imagenMuerte
let aud_muerte 
let aud_salto 
let aud_puntaje 
let generacion = 0

let layers = [3,3,1]
let num_genes = 13

let num_birds_gen = 100

var logo

class Bird{
    constructor(genoma){
        this.frames_1 = []
        let img = new Image()
        img.src = "sprites/yellowbird-downflap.png"
        this.frames_1.push(img)
        img = new Image()
        img.src = "sprites/yellowbird-midflap.png"
        this.frames_1.push(img)
        img = new Image()
        img.src = "sprites/yellowbird-upflap.png"
        this.frames_1.push(img)
        img = new Image()
        img.src = "sprites/yellowbird-midflap.png"
        this.frames_1.push(img)

        this.puntaje = 1
        
        this.pos = [40, 200]
        this.vel = [0, 0]
        this.acel = [0, 2500]
        
        this.jump_vel = [0,-600]
        
        this.rect = new Rect(this.pos, 32,27)
        
        this.image = this.frames_1[0];

        this.vida = 1

        if(genoma === undefined){
            genoma = new Genome(layers, num_genes, true, true)
    
            // Probabilidad de cambiar los parámetros como los sesgos y el peso de la conexión
            genoma.mutation_rate_parameter =  0.3/100;
            genoma.mutation_rate_connection = 0.1/100
            genoma.cruce = 70/100
        }
        this.cerebro = create_neural_network(genoma)
        this.cerebro.capas_ocultas = activation_functions.sigmoid
        this.cerebro.capa_salida = activation_functions.sigmoid
    }
    think(pipe){
        let input = [pipe.pos[1]-pipe.separacion/2,pipe.pos[1]+pipe.separacion/2, this.pos[1]]
        let output = this.cerebro.forward(input)
        if(output[0]>0.5){
            this.jump()
        }
    }
    jump(){
        if(this.acel[0]==0 && this.vida == 1){
            this.acel = [0, 4000]
        }
        aud_salto.currentTime = 0;
        aud_salto.play()
        this.vel = [...this.jump_vel]
    }
    updateFrames(){
        let tiempo = Date.now();
        this.image = this.frames_1[parseInt((tiempo/70)%4)]
    }
}

class Pipe{
    constructor(type, x, separacion, max_y){
        this.image1 = new Image()

        const self = this;
        
        if(type){
            this.image1.src = "sprites/pipe-red.png"
        }else{
            this.image1.src = "sprites/pipe-green.png"
        }

        this.image1.onload = function() {
            self.separacion = separacion
            self.max_y = max_y
            self.rect1 = new Rect([0,0], 47,320)
            self.rect2 = new Rect([0,0], 47,320)
            self.generar_pos(x)
        };
    }
    generar_pos(x){
        this.pos = [x, getRandomNumber(this.separacion/2,this.max_y-this.separacion/2)]
        
        if(this.rect1 !==undefined){
            this.rect1.pos = [this.pos[0],this.pos[1]-this.separacion/2-this.image1.height/2]
            this.rect2.pos = [this.pos[0],this.pos[1]+this.separacion/2+this.image1.height/2]
        }
    }
}

function generate_population(agentes, num_agentes){
    // Calcular la mitad del tamaño del arreglo original
    let mitad = Math.floor(agentes.length / 2);

    for (let i = mitad; i < agentes.length; i++) {
        const b = agentes[i];
        birds.push(new Bird(b.cerebro.genome))
    }


    while(birds.length<num_agentes){
        parents = parent_selection(agentes)

        children = parents[0].cerebro.genome.crossover(JSON.parse(JSON.stringify(parents[1].cerebro.genome)))
		children.forEach(son => {
			son.mutate() 
            birds.push(new Bird(son))
		});
    }
}

function setup(){

    canvas = document.getElementById('gameCanvas')
    ctx = canvas.getContext('2d')


    logo = new Image()
    logo.src = "logo.jpeg"

    
    birds = []
    for (let i = 0; i < num_birds_gen; i++) {
        genoma = new Genome(layers, num_genes, false, false)
        genoma.mutation_rate_parameter = 0.1
        genoma.cruce = 60/100
        birds.push(new Bird(genoma))
    }
    
    for (let i = 0; i < 3; i++) {
        let image = new Image()
        image.src = "sprites/background-day.png"
        fondo.push(image)
    }
    for (let i = 0; i < 5; i++) {
        let image = new Image()
        image.src = "sprites/base.png"
        suelo.push(image)
    }
    
    for (let i = 0; i < 6; i++) {
        tubos.push(new Pipe(Math.random()>0.8, canvas.width/2+i*distancia_regeneracion, 170, canvas.height-suelo[0].height))
    }

    imagenMuerte = new Image()
    imagenMuerte.src = "sprites/gameover.png"
    
    
    rectSuelo = new Rect([0,canvas.height-suelo[0].height/2], 100, 110)

    aud_muerte = new Audio('audio/hit.ogg');
    aud_salto = new Audio('audio/wing.ogg');
    aud_puntaje = new Audio('audio/point.ogg');

    velocidad_juego = -200

    birds_Eliminados = []
}

function update(dt){
    if (velocidad_juego>-600){
        velocidad_juego -=0.2
    }
    
    if(birds.length>0){
        draw_neural_network(birds[0].cerebro)
    } 
    
    if(getKeyDown(keyCodes.Space)){
        
        // birds.forEach(bird => {
        //     bird.jump()
        // });
        console.log(JSON.stringify(birds[0].cerebro.biases))
        console.log(JSON.stringify(birds[0].cerebro.weights))
    }
    
    fondo.forEach((img, index) => {
        drawImage(ctx,img, img.width*index+img.width/2,img.height/2,0)
    });
    
    suelo_x+=velocidad_juego*dt
    if(suelo_x<=-suelo[0].width){
        suelo_x=0
    }

    tubos.forEach(tubo => {
        tubo.pos[0]+=velocidad_juego*dt
        
        tubo.rect1.pos[0]+=velocidad_juego*dt
        tubo.rect2.pos[0]+=velocidad_juego*dt
        drawImage(ctx, tubo.image1, tubo.pos[0],tubo.pos[1]+tubo.separacion/2+tubo.image1.height/2, 0)
        drawImage(ctx, tubo.image1, tubo.pos[0],tubo.pos[1]-tubo.separacion/2-tubo.image1.height/2, 180*radians)
        
        
        birds.forEach(bird => {
            if(tubo.rect1.isCollidingWith(bird.rect) && bird.vida>0){
                bird.vida = 0
                aud_muerte.play()
            }
            if(tubo.rect2.isCollidingWith(bird.rect) && bird.vida>0){
                bird.vida = 0
                aud_muerte.play()
            }
        });
    });
    suelo.forEach((img, index) => {
        drawImage(ctx,img, img.width*index+img.width/2+suelo_x,canvas.height-img.height/2,0)
    });
    if(tubos[0].pos[0]<-tubos[0].image1.width){
        const firstElement = tubos.shift();
        firstElement.generar_pos(tubos[tubos.length-1].pos[0]+distancia_regeneracion)
        tubos.push(firstElement);
        puntaje+=1
        
        aud_puntaje.currentTime = 0;
        aud_puntaje.play()
    }

    birds.forEach(bird => {
        bird.puntaje+=0.01
    });
    
    
    birds.forEach((bird, i) => {
        if(i<20){
            drawImage(ctx, bird.image, bird.pos[0], bird.pos[1], (bird.vel[1]/reduccion_giro)*radians)
        }
        // if (i==0)
        bird.vel[1]+=bird.acel[1]*dt
        bird.pos[1]+=bird.vel[1]*dt
        bird.updateFrames()
        if(rectSuelo.isCollidingWith(bird.rect) && bird.vida>0){
            bird.vida = 0
        }

        if(bird.pos[1]<-20 || bird.pos[1]>canvas.height){
            bird.vida=0
        }
        
        if(bird.vida == 0){
            const removedElement = birds.splice(i, 1)[0];
            birds_Eliminados.push(removedElement)
        }else{
            bird.think(tubos[tubos[0].pos[0]>30?0:1])
        }
    });
    if(birds.length==0){
        velocidad_juego = -200
        
        generate_population(birds_Eliminados, num_birds_gen)
        birds_Eliminados=[]
        tubos.forEach((tubo, i) => {
            tubo.generar_pos(canvas.width/2+i*distancia_regeneracion)
        });
        puntaje = 0
        generacion +=1
    }

    drawImageWithOpacity(ctx, logo, canvas.width/2, canvas.height/2-70, 0,0.4, 100,100)

    drawNumber(ctx, puntaje, canvas.width/2, 50)
    drawNumber(ctx, birds.length, canvas.width/10, 20)
    drawNumber(ctx, generacion, canvas.width/10*9, 20)

    // drawRect(ctx, tubos[0].rect1, "red")
}
function gameLoop() {
    let lastTime = 0;
    let delay = 1000;
    let startTime = null;
    

    function loop(timestamp) {
        if (!startTime) {
            startTime = timestamp;
        }

        const elapsed = timestamp - startTime;
        

        if (elapsed < delay) {
            lastTime = timestamp;
            requestAnimationFrame(loop);
            return;
        }

        const dt = (timestamp - lastTime) / 1000; 
        lastTime = timestamp;
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        
        update(dt);
        requestAnimationFrame(loop);
    }
    
    requestAnimationFrame(loop);
}

setup()
gameLoop();
