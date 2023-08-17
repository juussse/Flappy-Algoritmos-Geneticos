function drawRectangle(context, x, y, width, height, color) {
    context.fillStyle = color;
    context.fillRect(x, y, width, height);
}
function drawRect(context, rect, color){
    drawRectangle(ctx, rect.pos[0]-rect.width/2,rect.pos[1]-rect.height/2, rect.width, rect.height, color)
}
function drawTriangle(context, x1, y1, x2, y2, x3, y3, color) {
    context.beginPath();
    context.moveTo(x1, y1);
    context.lineTo(x2, y2);
    context.lineTo(x3, y3);
    context.closePath();
    context.fillStyle = color;
    context.fill();
}
function drawLine(context, x1, y1, x2, y2, color, lineWidth) {
    context.beginPath();
    context.moveTo(x1, y1);
    context.lineTo(x2, y2);
    context.strokeStyle = color;
    context.lineWidth = lineWidth;
    context.stroke();
}
function drawCircle(context, x, y, radius, color, borderWidth) {
    context.beginPath();
    context.arc(x, y, radius, 0, 2 * Math.PI);
    
    if (borderWidth > 0) {
        context.strokeStyle = color;
        context.lineWidth = borderWidth;
        context.stroke();
    } else {
        context.fillStyle = color;
        context.fill();
    }
    
    context.closePath();
}
function drawImage(context, image, x, y, angle, nWidth, nHeight) {
    const width = nWidth===undefined?image.width:nWidth;
    const height = nHeight===undefined?image.height:nHeight;
    
    // Calcular las coordenadas superiores izquierdas del rectángulo que contenga la imagen centrada en (x, y)
    const topLeftX = x - width / 2;
    const topLeftY = y - height / 2;
    
    // Rotar el contexto según el ángulo especificado
    context.save();
    context.translate(x, y);
    context.rotate(angle);
    
    // Dibujar la imagen centrada
    context.drawImage(image, -width / 2, -height / 2, width, height);
    
    // Restaurar el contexto
    context.restore();
}
function drawImageWithOpacity(context, image, x, y, angle, opacity, nWidth, nHeight) {
    const width = nWidth === undefined ? image.width : nWidth;
    const height = nHeight === undefined ? image.height : nHeight;

    // Calcular las coordenadas superiores izquierdas del rectángulo que contenga la imagen centrada en (x, y)
    const topLeftX = x - width / 2;
    const topLeftY = y - height / 2;

    // Rotar el contexto según el ángulo especificado
    context.save();
    context.translate(x, y);
    context.rotate(angle);

    // Establecer la opacidad
    context.globalAlpha = opacity;

    // Dibujar la imagen centrada con la opacidad especificada
    context.drawImage(image, -width / 2, -height / 2, width, height);

    // Restaurar el contexto y la opacidad
    context.restore();
    context.globalAlpha = 1.0;
}


// Función para dibujar un puntaje en el lienzo utilizando imágenes individuales para cada dígito
function drawNumber(ctx, score, x, y) {
    // Convierte el puntaje a una cadena de caracteres
    var scoreString = score.toString();
    
    // Tamaño de los dígitos en píxeles
    var digitWidth = 20;
    var digitHeight = 30;
    
    // Ciclo para dibujar cada dígito del puntaje
    for (var i = 0; i < scoreString.length; i++) {
        // Obtén el dígito actual
        var digit = parseInt(scoreString.charAt(i));
        
        // Crea una nueva imagen para el dígito actual
        var digitImage = new Image();
        digitImage.src = "sprites/"+ digit + '.png'; // Nombre de archivo correspondiente al dígito
        
        // Dibuja el dígito en el lienzo
        ctx.drawImage(digitImage, x-(scoreString.length*digitWidth)/2, y, digitWidth, digitHeight);
        
        // Aumenta la posición x para el siguiente dígito
        x += digitWidth;
    }
}

function drawText(context, text, x, y, fontSize, fontColor, fontFamily, textAlign) {
    context.font = fontSize + 'px ' + fontFamily;
    context.fillStyle = fontColor;
    context.textAlign = textAlign;
    context.fillText(text, x, y);
}

const canvasRed = document.getElementById("info")
const ctxR = canvasRed.getContext('2d')
function draw_neural_network(neural_n){
    
    ctxR.clearRect(0,0, canvasRed.width, canvasRed.height)
    let max_num_neu = 0
    for(let i=0;i<neural_n.layer_sizes.length; i++){
        if (max_num_neu<neural_n.layer_sizes[i]){
            max_num_neu = neural_n.layer_sizes[i]
        }
    } 
    
    let sepX = canvasRed.width/(neural_n.layer_sizes.length+1)
    let sepY = canvasRed.height/(max_num_neu+1)
    
    radius = Math.min(30, Math.min(sepX, sepY)/2-(Math.min(sepX, sepY)*0.1))
    
    
    
    for(let i=0;i<neural_n.layer_sizes.length; i++){
        n_neu = neural_n.layer_sizes[i]
        for(let j=0; j<n_neu;j++){
            cN1 = [sepX*i+sepX, sepY*j+canvasRed.height/2-(sepY*(neural_n.layer_sizes[i]-1))/2]
            
            drawCircle(ctxR, cN1[0],cN1[1], radius, "red", 1)
            
            if (i>0){
                b = neural_n.biases[(i-1)][j]
                drawText(ctxR, b.toFixed(2), cN1[0], cN1[1]+3, 10,"black", "sans-serif", "center")
            }
            
            if (i<neural_n.layer_sizes.length-1){
                for(let k=0; k<neural_n.layer_sizes[i+1]; k++){
                    cN2 = [sepX*(i+1)+sepX, sepY*k+canvasRed.height/2-(sepY*(neural_n.layer_sizes[i+1]-1))/2]
                    w = neural_n.weights[i][k][j]
                    // print("PESOS "+w)
                    if (Math.abs(w)>0){
                        drawLine(ctxR, cN1[0]+radius, cN1[1], cN2[0]-radius, cN2[1], "red", 1, Math.abs(w)*200)
                    }
                }
            }
        }
    }
}

function dibujarCirculoDifuminado(ctx, x, y, radio, colorCentro, colorBorde, colorBordeExt) {
    const gradiente = ctx.createRadialGradient(x, y, 0, x, y, radio);
    gradiente.addColorStop(0, colorCentro);    // Color en el centro
    gradiente.addColorStop(1, colorBorde);     // Color en el borde
    gradiente.addColorStop(1, colorBordeExt);     // Color en el borde
    
    ctx.beginPath();
    ctx.arc(x, y, radio, 0, 2 * Math.PI);
    ctx.fillStyle = gradiente;
    ctx.fill();
    ctx.closePath();
}


keys = {};
keysDown = {};

function keyDownHandler(event) {
    const keyCode = event.keyCode || event.which;
    keys[keyCode] = true;
    keysDown[keyCode] = true
}

function getKeyDown(code){
    if (keysDown[code]){
        keysDown[code] = false
        return true
    }
    return false
}
function keyUpHandler(event) {
    const keyCode = event.keyCode || event.which;
    keys[keyCode] = false;
}


window.addEventListener('keydown', keyDownHandler);
window.addEventListener('keyup', keyUpHandler);