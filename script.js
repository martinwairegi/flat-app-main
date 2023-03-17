let RADIUS = 40;        // state radius
let CHEVRON = RADIUS/4; // length of transition chevron
let SELECTAREA = 10;    // padding either side of transitions for easier selection
let FONTSIZE = 16;      // font size for labels
let EPSILON = String.fromCharCode(949); // epsilon symbol
let SIGMA = ['a','b', 'c'];  // fsm alphabet
let STATEFILL = "#fcfcfc" // fill colour of states
let BLACK = "#000000"   // black hex code
let RED = "#ff0000"     // red hex code
const nodes = [];       // array of states
var edges = [];         // array of transitions
var sid = 0;            // unique state ID
var tid = 0;            // unique transition ID
var highSid = -1;       // ID of highlighted state
var highTid = -1;       // ID of highlighted transition
var startSid = -1;      // ID of start state
var startTid = -1;      // ID of start transition

class Regex {

    constructor() {
        this.generate();
    }

    /**
     * Generate new regular expression and corresponding NFA
     */
    generate() {
        this.postfix = "";
        this.regex = this.#kleene(7, 0.5, 0.2, 0.1);
        this.nfa = this.#regexToNfa(this.postfix);
    }

    /**
     * Convert regular expression into equivalent NFA using Thompson's construction
     * @param {String} regex A regular expression in postfix
     * @returns {
     *  Array<{
     *      table: 
     *          Array<{
     *              stateID: Number,
     *              symbol: String,
     *              stateIDs: Array<Number>
     *          }>,
     *      start: Number,
     *      end: Array<Number>
     *  }>
     * } State transition table for NFA accepting regex
     *  (alongside start state and accept states)
     */
    #regexToNfa(regex) {
        const nfa = []; // State transition table
        const s = [];   // Stack of pairs of states to next consider
        var start = 0;  // ID of start state
        var end = 1;    // ID of accept state
        var count = 0;  // Counter for state IDs
        var c1 = 0;     // ID of a state to add to NFA
        var c2 = 0;     // ID of another state to add to NFA

        // Iterate through each character in the postfix regex
        for (var i=0; i<regex.length; i++) {
            if (regex[i] == '*') { // Kleene star operator
                // Pop last pair of states from stack (sub-NFA)
                var top = s.pop();
                var r1 = top[0]; // start of sub-NFA
                var r2 = top[1]; // end of sub-NFA
                // Set IDs of two new states
                c1 = count++;
                c2 = count++;
                // Push new states onto stack
                s.push([c1, c2]);
                // Add new states to NFA
                nfa.push({});
                nfa.push({});
                for (var char of SIGMA) {
                    nfa[c1][char] = [];
                    nfa[c2][char] = [];
                }
                nfa[c1][EPSILON] = [];
                nfa[c2][EPSILON] = []
                // Loop back to start of sub-NFA or continue
                nfa[r2][EPSILON].push(r1, c2);
                // Go to start of sub-NFA or skip
                nfa[c1][EPSILON].push(r1, c2);
                // Set new start and end states if necessary
                if (start == r1) {
                    start = c1;
                }
                if (end == r2) {
                    end = c2;
                }
            } else if (regex[i] == '.') { // Concatenation operator
                // Pop last two pairs of states from stack (two sub-NFAs)
                var top1 = s.pop();
                var top2 = s.pop();
                var r11 = top1[0];
                var r12 = top1[1];
                var r21 = top2[0];
                var r22 = top2[1];
                // Push 'start' of second pair and 'end' of first pair onto stack
                s.push([r21, r12]);
                // Connect first sub-NFA to second with epsilon transition
                nfa[r22][EPSILON].push(r11);
                // Set new start and end states if necessary
                if (start == r11) {
                    start = r21;
                }
                if (end == r22) {
                    end = r12;
                }
            } else if (regex[i] == '+') { // Or operator
                // Set IDs of two new states and add to NFA
                c1 = count++;
                c2 = count++;
                nfa.push({});
                nfa.push({});
                for (var char of SIGMA) {
                    nfa[c1][char] = [];
                    nfa[c2][char] = [];
                }
                nfa[c1][EPSILON] = [];
                nfa[c2][EPSILON] = []
                // Pop last two pairs of states from stack (two sub-NFAs)
                var top1 = s.pop();
                var top2 = s.pop();
                var r11 = top1[0];
                var r12 = top1[1];
                var r21 = top2[0];
                var r22 = top2[1];
                // Push new states to stack
                s.push([c1,c2]);
                // Traverse to second sub-NFA or first sub-NFA
                nfa[c1][EPSILON].push(r21, r11);
                // Continue from end of first sub-NFA
                nfa[r12][EPSILON].push(c2);
                // Continue from end of second sub-NFA
                nfa[r22][EPSILON].push(c2);
                // Set new start and end states if necessary
                if (start == r11 || start == r21) {
                    start = c1;
                }
                if (end == r22 || end == r12) {
                    end = c2;
                }
            } else { // symbol read
                // Set IDs of two new states and add to NFA
                c1 = count++;
                c2 = count++;
                nfa.push({});
                nfa.push({});
                for (var char of SIGMA) {
                    nfa[c1][char] = [];
                    nfa[c2][char] = [];
                }
                nfa[c1][EPSILON] = [];
                nfa[c2][EPSILON] = []
                // Push new states onto stack
                s.push([c1,c2]);
                // Connect the first state to the second via the symbol
                nfa[c1][regex[i]].push(c2);
            }
        }

        return {
            "table" : nfa,
            "start" : start,
            "end" : end
        }
    }

    /**
     * Closes expr in the Kleene star (*) operator with probability probKleene
     * @param {Number} n Number of terms in regex
     * @param {Number} probOr Probability Or (+) operator used in favour of Concatenation (.) operator
     * @param {Number} probKleene Probability Kleene star (*) used
     * @param {Number} probEmpty Probability epsilon character used
     * @returns {String} Regular expression
     */
    #kleene(n, probOr, probKleene, probEmpty) {
        // Generate expression
        var expr = this.#expression(n, probOr, probKleene, probEmpty);
        // Apply Kleene star operator with probability probKleene
        if (Math.random() <= probKleene) {
            if (expr.length > 1) {
                expr = "(" + expr + ")*";
            } else {
                expr = expr + "*";
            }
            this.postfix += "*";
        }
        return expr;
    }

    /**
     * Constructs a regular expression with operators and symbols included probabilistically
     * @param {Number} n Number of terms in regex
     * @param {Number} probOr Probability Or (+) operator used in favour of Concatenation (.) operator
     * @param {Number} probKleene Probability Kleene star (*) used
     * @param {Number} probEmpty Probability epsilon character used 
     * @returns {String} Regular expression
     */
    #expression(n, probOr, probKleene, probEmpty) {
        if (n < 2) {
            // Randomly select symbol from sigma
            var symbol = SIGMA[Math.floor(Math.random() * SIGMA.length)];
            this.postfix += symbol;
            return symbol;
        } else if (Math.random() <= probEmpty) { // use epsilon with probability probEmpty
            this.postfix += EPSILON;
            // Generate smaller sub-expression
            var after = this.#kleene(n-1, probOr, probKleene, probEmpty);
            this.postfix += "+";
            return "(" + EPSILON + " + " + after + ")";
        }

        var beforeSize = Math.floor(n/2);

        // Generate two sub-expressions
        var before = this.#kleene(beforeSize, probOr, probKleene, probEmpty);
        var after = this.#kleene(n-beforeSize, probOr, probKleene, probEmpty);

        // Apply Or operator between the two with probability probOr
        if (Math.random() <= probOr) {
            this.postfix += "+";
            return "(" + before + " + " + after + ")";
        }

        // Apply Concatenation operator between the two with probability 1-probOr
        this.postfix += ".";
        return before + after;
    }

}

class Edge {

    constructor(id, fromNode, toNode) {
        this.id = id;
        this.fromNode = fromNode;
        this.toNode = toNode;
        this.label = "";

        // Set if self loop
        this.x = null;
        this.y = null;
        this.radius = null;

        // Set if non self loop
        this.angle = null;

        // Set if curved
        this.curved = false;
    }

    /**
     * Draws edge to canvas
     * @param {CanvasRenderingContext2D} ctx 2D rendering context for drawing surface of FSM canvas
     */
    draw(ctx) {

        ctx.strokeStyle = BLACK;
        ctx.fillStyle = BLACK;

        // Colour edge red if highlighted
        if (this.id == highTid) {
            ctx.strokeStyle = RED;
            ctx.fillStyle = RED;
        }

        ctx.beginPath();

        if (this.fromNode == this.toNode) { // self loop
            this.angle = 5*Math.PI/16;
            var dx = Math.cos(this.angle)*RADIUS;
            var dy = Math.sin(this.angle)*RADIUS;
            var xn = this.fromNode.x;
            var yn = this.fromNode.y;

            // Start of arc
            var x1 = xn-dx;
            var y1 = yn-dy;
            // End of arc
            var x2 = xn+dx;
            var y2 = yn-dy;
            // Arc turning point
            var x3 = xn;
            var y3 = yn-1.7*RADIUS;

            // Find circle equation from three points (above)
            var circle = circleFromPoints(x1, y1, x2, y2, x3, y3);

            this.x = circle.x; // x centre
            this.y = circle.y // y centre
            this.radius = circle.radius;

            // Angle between arc centre and end of arc
            var alpha = Math.atan2(y2-this.y, x2-this.x); 

            ctx.beginPath();
            ctx.arc(this.x, this.y, this.radius, Math.PI-alpha, alpha); // arc is drawn outside of node area
            ctx.stroke();

            // Draw chevron at end of arc
            ctx.beginPath();
            ctx.moveTo(x2, y2);
            ctx.lineTo(x2+CHEVRON*Math.cos(this.angle-Math.PI/10), y2-CHEVRON*Math.sin(this.angle-Math.PI/10));
            ctx.lineTo(x2-CHEVRON*Math.cos(this.angle+Math.PI/10), y2-CHEVRON*Math.sin(this.angle+Math.PI/10));
            ctx.closePath();
            ctx.stroke();
            ctx.fill();

            ctx.strokeStyle = BLACK; // revert colour to black

            ctx.fillStyle = STATEFILL;
                
            var width = ctx.measureText(this.label).width;

            ctx.fillRect(x3-width/2, y3-4-FONTSIZE+2, width, FONTSIZE+2);

            ctx.fillStyle = BLACK;

            ctx.beginPath();
            ctx.fillText(this.label, x3, y3-4);
            ctx.stroke();

            ctx.fillStyle = STATEFILL
        } else if (this.curved) { // curved edge between nodes
            var x1 = this.fromNode.x;
            var y1 = this.fromNode.y;

            var x2 = this.toNode.x;
            var y2 = this.toNode.y;

            var dx = x1-x2;
            var dy = y1-y2;
            
            this.angle = Math.atan2(dy, dx);

            var x3 = 0.5*(x1+x2) + 2*SELECTAREA*Math.cos(this.angle - Math.PI/2);
            var y3 = 0.5*(y1+y2) + 2*SELECTAREA*Math.sin(this.angle - Math.PI/2);

            // create circle using three points
            var circle = circleFromPoints(x1, y1, x2, y2, x3, y3);

            var xc = circle.x;
            var yc = circle.y;

            // only draw section between nodes
            var startAngle = Math.atan2(y2-yc, x2-xc);
            var endAngle = Math.atan2(y1-yc, x1-xc);

            ctx.beginPath();
            ctx.arc(xc, yc, circle.radius, startAngle, endAngle);
            ctx.stroke();

            // get coords of arc intersection with 'to' node
            var alpha = Math.acos(RADIUS/(2*circle.radius)) - startAngle + Math.PI;

            var xi = x2 + RADIUS*Math.cos(alpha);
            var yi = y2 - RADIUS*Math.sin(alpha);

            var beta = Math.atan2(yi-y2,xi-x2);
            
            // dynamically draw chevron
            ctx.beginPath();
            ctx.moveTo(xi, yi);
            ctx.lineTo(xi+CHEVRON*Math.cos(beta-Math.PI/5), yi+CHEVRON*Math.sin(beta-Math.PI/5));
            ctx.lineTo(xi+CHEVRON*Math.cos(beta+Math.PI/5), yi+CHEVRON*Math.sin(beta+Math.PI/5));
            ctx.closePath();
            ctx.stroke();
            ctx.fill();

            ctx.strokeStyle = BLACK; // revert colour to black

            // draw the label at the third point that was created
            ctx.fillStyle = STATEFILL;
                
            var width = ctx.measureText(this.label).width;

            ctx.fillRect(x3-width/2, y3-FONTSIZE+2, width, FONTSIZE+2);

            ctx.fillStyle = BLACK;

            ctx.beginPath();
            ctx.fillText(this.label, x3, y3);
            ctx.stroke();

            ctx.fillStyle = STATEFILL;
        } else {
            if (this.id == startTid) { // start edge
                var toX = this.toNode.x-RADIUS;
                var toY = this.toNode.y;
                var fromX = toX-RADIUS;
                var fromY = toY;
                var dx = RADIUS;
                var dy = 0;
                this.angle = Math.atan2(dy, dx);
            } else { // edge between nodes
                var toX = this.toNode.x;
                var toY = this.toNode.y;
                var fromX = this.fromNode.x;
                var fromY = this.fromNode.y;

                // Calculates line angle between centres of each node
                var dx = toX-fromX;
                var dy = toY-fromY;
                this.angle = Math.atan2(dy, dx);

                // 'Remove' portion of edge contained within nodes
                fromX += Math.cos(this.angle)*RADIUS;
                fromY += Math.sin(this.angle)*RADIUS;
                toX -= Math.cos(this.angle)*RADIUS;
                toY -= Math.sin(this.angle)*RADIUS;
            }

            // Draw connecting line
            ctx.beginPath();
            ctx.moveTo(fromX, fromY);
            ctx.lineTo(toX, toY);
            ctx.stroke();

            // Draw chevron at end of edge
            ctx.beginPath();
            ctx.moveTo(toX, toY);
            ctx.lineTo(toX-CHEVRON*Math.cos(this.angle - Math.PI/6), toY-CHEVRON*Math.sin(this.angle - Math.PI/6));
            ctx.lineTo(toX-CHEVRON*Math.cos(this.angle + Math.PI/6), toY-CHEVRON*Math.sin(this.angle + Math.PI/6));
            ctx.closePath();
            ctx.stroke();
            ctx.fill();

            ctx.strokeStyle = BLACK; // revert colour to black
            ctx.fillStyle = STATEFILL;

            if (this.fromNode != null) {

                var width = ctx.measureText(this.label).width;

                var x = (this.fromNode.x + this.toNode.x) / 2;
                var y = (this.fromNode.y + this.toNode.y) / 2;

                ctx.fillRect(x-width/2, y-FONTSIZE+2, width, FONTSIZE+2);

                ctx.fillStyle = BLACK;

                ctx.beginPath();
                ctx.fillText(this.label, x, y);
                ctx.stroke();

                ctx.fillStyle = STATEFILL;
            }
        }
    }
}


class Node {

    constructor(id, x, y) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.label = "";
        this.accept = false;
        this.dragging = false;
    }

    /**
     * Draws node to canvas
     * @param {CanvasRenderingContext2D} ctx 2D rendering context for drawing surface of FSM canvas
     */
    draw(ctx) {
        // Colour state red if highlighted
        if (this.id == highSid) {
            ctx.strokeStyle = RED;
        }

        // Draw state
        ctx.beginPath();
        ctx.arc(this.x, this.y, RADIUS, 0, 2*Math.PI);
        ctx.fill();
        ctx.stroke();

        // Draw smaller circle inside to denote accept state
        if (this.accept) {
            ctx.beginPath();
            ctx.arc(this.x, this.y, RADIUS-8, 0, 2*Math.PI);
            ctx.fill();
            ctx.stroke();
        }

        ctx.strokeStyle = BLACK; // revert colour to black

        ctx.fillStyle = BLACK;
        ctx.beginPath();
        ctx.fillText(this.label, this.x, this.y+5);
        ctx.stroke();

        ctx.fillStyle = STATEFILL;
    }
}

/**
 * Determines if DFAs given by 'user' and 'regex' accept the same language
 * @param {
 *  Array<{
 *      dfa: 
 *          Array<{
 *              stateIdFrom: Number,
 *              symbol: String,
 *              stateIdTo: Number
 *          }>,
 *      start: Number,
 *      accept: Array<Number>
 *  }>
 * } user State transition table, start state, and accept states defined for the user's DFA
 * @param {
 *  Array<{
 *      dfa: 
 *          Array<{
 *              stateIdFrom: Number,
 *              symbol: String,
 *              stateIdTo: Number
 *          }>,
 *      start: Number,
 *      accept: Array<Number>
 *  }>
 * } regex State transition table, start state, and accept states defined for the regex's DFA
 * @returns {Boolean} True iff user DFA and regex DFA accept the same language
 */
function isomorphic(user, regex) {
    // Array of all accept states (state IDs are unique across both DFAs)
    const accept = user.accept.concat(regex.accept);
    // Each state has a set represented by a tree: 'parent' represents the root
    //  of the tree for a given state
    const parent = {};
    // Each state also has a rank denoting their height in their tree
    const rank = {};
    // Stack containing pairs of states
    const pairStack = [];
    // User states
    const userStates = [];
    for (var k of Object.keys(user.dfa)) {
        userStates.push(parseInt(k));
    }

    // Make sets for each state (represented as trees)
    for (var s of Object.keys(user.dfa)) {
        makeSet(parseInt(s), parent, rank);
    }
    for (var s of Object.keys(regex.dfa)) {
        makeSet(parseInt(s), parent, rank);
    }

    // Assume DFAs are equal to begin with
    var equal = true;

    // Calculate union of start states
    equal = unionCheck(user.start, regex.start, parent, rank, accept);

    // Push start states of both DFAs to stack
    pairStack.push([user.start, regex.start]);

    // While the stack is nonempty and condition of equivalence has not yet been violated
    while (pairStack.length > 0 && equal) {
        // Pop next pair of states
        pair = pairStack.pop();
        // For each symbol in the alphabet
        for (var c of SIGMA) {
            // Take transition via 'c' for each state and determine which set they belong to
            //  (i.e., the root of the tree they're in)
            var r1 = 0;
            var r2 = 0;
            if (userStates.includes(pair[0])) {
                r1 = findSet(user.dfa[pair[0]][c], parent);
            } else {
                r1 = findSet(regex.dfa[pair[0]][c], parent);
            }
            if (userStates.includes(pair[1])) {
                r2 = findSet(user.dfa[pair[1]][c], parent);
            } else {
                r2 = findSet(regex.dfa[pair[1]][c], parent);
            }
            // If they belong to different sets
            if (r1 != r2) {
                // Take the union of the sets
                equal = unionCheck(r1, r2, parent, rank, accept);
                // Push the traversed to-states onto the stack
                pairStack.push([r1, r2]);
            }
        }
    }

    return equal;
}

/**
 * Create new tree rooted at state ID x (represents a set for x)
 * @param {Number} x State ID
 * @param {Array<Number>} parent Array of root nodes for each tree set
 * @param {Array<Number>} rank Array of rank for each node
 */
function makeSet(x, parent, rank) {
    parent[x] = x;
    rank[x] = 0;
}

/**
 * Calculates the union of sets x and y and whether the user and regex
 *  DFAs are isomorphic up to this point
 * @param {Number} x State ID
 * @param {Number} y State ID
 * @param {Array<Number>} parent Array of root nodes for each tree set
 * @param {Array<Number>} rank Array of rank for each node
 * @param {Boolean} accept Flag to see if user and regex DFAs are isomorphic
 * @returns {Boolean} True iff the union of sets x and y contain either
 *  only accepting states, or only non-accepting states
 */
function unionCheck(x, y, parent, rank, accept) {
    // Get root nodes for each state
    var a = findSet(x, parent);
    var b = findSet(y, parent);
    // Return false if the union of sets contains both a non-accepting
    //  and accepting state
    if (accept.includes(a)) {
        if (!accept.includes(b)) {
            return false;
        }
    } else {
        if (accept.includes(b)) {
            return false;
        }
    }
    // If the union contains either only accepting states, or only
    //  non-accepting states, link them and return true
    link(a, b, parent, rank);
    return true;
}

/**
 * Sets the root node of the lower ranked state as the higher ranked state
 * @param {Number} x State ID
 * @param {Number} y State ID
 * @param {Array<Number>} parent Array of root nodes for each tree set
 * @param {Array<Number>} rank Array of rank for each node
 */
function link(x, y, parent, rank) {
    // Set higher ranked state as root node of lower ranked state
    if (rank[x] > rank[y]) {
        parent[y] = x;
    } else {
        parent[x] = y;
        if (rank[x] == rank[y]) {
            rank[y] += 1;
        }
    }
}

/**
 * Finds the root node for state x
 * @param {Number} x State ID
 * @param {Array<Number>} parent Array of root nodes for each tree set
 * @returns {Number} Root node for state x
 */
function findSet(x, parent) {
    if (x != parent[x]) {
        // Recursively find parent of x until root is reached
        parent[x] = findSet(parent[x], parent);
    }
    return parent[x];
}

/**
 * Convert NFA to a DFA using subset construction
 * @param {
 *  Array<{
 *      stateID: Number,
 *      symbol: String,
 *      stateIDs: Array<Number>
 *  }>
 * } nfa State transition table
 * @param {Number} start NFA start state
 * @param {Array<Number>} final NFA accept states
 * @param {Number} dfaId ID of initial DFA state
 * @returns {
 *  Array<{
 *      dfa: 
 *          Array<{
 *              stateIdFrom: Number,
 *              symbol: String,
 *              stateIdTo: Number
 *          }>,
 *      start: Number,
 *      accept: Array<Number>
 *  }>
 * } State transition table, start state, and accept states defined for the resulting DFA
 */
function subsetConstruct(nfa, start, final, dfaId) {
    // Accept states of DFA
    const accept = [];

    // DFA start state
    var begin = dfaId;

    // State transition table
    const dfa = {};
    // Mapping of state subsets in NFA to state IDs in DFA
    const dfaIds = {};

    // Contains E-CLOSE for each state
    const nodeClosure = [];
    for (const [n, t] of Object.entries(nfa)) {
        nodeClosure[n] = [];
    }

    // First state of DFA by performing E-CLOSE on start state of NFA
    var firstState = eClose([start], nodeClosure, nfa);
    // Add new state to DFA
    dfa[dfaId] = {};
    // Map subset given by E-CLOSE to corresponding DFA state ID
    dfaIds[firstState] = dfaId++;
    // DFA state is accepting if corresponding state subset contains accepting state
    for (var n of firstState) {
        if (final.includes(n)) {
            accept.push(dfaIds[firstState]);
            break;
        }
    }

    // Initialise a queue of state subsets to process
    const nodeQueue = [firstState];
    
    while (nodeQueue.length > 0) {
        // Dequeue next state subset
        var currentState = nodeQueue.shift();
        // For each symbol of the alphabet
        for (var s of SIGMA) {
            // State subset calculated by reading 's' from each node in current subset
            //  and applying E-CLOSE to the result
            var subset = nodeSubset(currentState, s, nodeClosure, nfa).sort();
            // If result is a new subset of states
            if (!(subset in dfaIds)) {
                // Add new state to DFA
                dfa[dfaId] = {};
                // Map subset to corresponding DFA state ID
                dfaIds[subset] = dfaId++;
                // DFA state is accepting if corresponding state subset contains accepting state
                for (var n of subset) {
                    if (final.includes(n)) {
                        accept.push(dfaIds[subset]);
                        break;
                    }
                }
                // Queue new subset for processing
                nodeQueue.push(subset);
            }
            // Add transition between the state subsets
            dfa[dfaIds[currentState]][s] = dfaIds[subset];
        }
    }

    return {
        dfa : dfa,
        start : begin,
        accept : accept
    }
}

/**
 * Calculates subset of states reached by reading 'symbol' from each state
 *  in 'states' then applying E-CLOSE to it
 * @param {Array<Number>} states Subset of states in NFA
 * @param {String} symbol State transition symbol
 * @param {Array<{stateID: Array<Number>}>} nodeClosure Closure set for each state in NFA
 * @param {
 *  Array<{
 *      stateID: Number,
 *      symbol: String,
 *      stateIDs: Array<Number>
 *  }>
 * } nfa State transition table
 * @returns {Array<Number>} Subset of NFA states that constitute a state in the DFA
 */
function nodeSubset(states, symbol, nodeClosure, nfa) {
    var subset = new Set();
    // For each state in subset of states in NFA
    for (var s of states) {
        // For each transition outgoing from state
        for (var t in nfa[s]) {
            // If transition symbols match
            if (t == symbol) {
                // Add each state reached by reading symbol
                for (var n of nfa[s][t]) {
                    subset.add(n);
                }
            }
        }
    }
    var nodeIds  = [];
    for (var n of subset.values()) {
        nodeIds.push(n);
    }
    // Return E-CLOSE of result
    return eClose(nodeIds, nodeClosure, nfa);
}

/**
 * Calculates E-CLOSE of 'states'
 * @param {Array<Number>} states Subset of states in NFA
 * @param {Array<{stateID: Array<Number>}>} nodeClosure Closure set for each state in NFA
 * @param {
 *  Array<{
 *      stateID: Number,
 *      symbol: String,
 *      stateIDs: Array<Number>
 *  }>
 * } nfa State transition table
 * @returns {Array<Number>} Subset of states as a result of applying E-CLOSE to 'states'
 */
function eClose(states, nodeClosure, nfa) {
    var closed = new Set();
    // For each state in subset of states
    for (var n of states) {
        // If closure for given state not yet calculated
        if (nodeClosure[n].length == 0) {
            var nClosed = new Set();
            // Calculate E-CLOSE for state
            var eStates = close(n, nodeClosure, nClosed, nfa);
            for (var q of eStates) {
                nodeClosure[n].push(q);
            }
        }
        // Add each state in closure of given state
        for (var q of nodeClosure[n]) {
            closed.add(q);
        }
    }
    const values = [];
    for (var v of closed.values()) {
        values.push(v);
    }
    return values;
}

/**
 * Calculates E-CLOSE for a specific state k
 * @param {Number} k State ID
 * @param {Array<{stateID: Array<Number>}>} nodeClosure Closure set for each state in NFA 
 * @param {Set<Number>} nClosed E-CLOSE set
 * @param {
 *  Array<{
 *      stateID: Number,
 *      symbol: String,
 *      stateIDs: Array<Number>
 *  }>
 * } nfa State transition table
 * @returns {IterableIterator<Number>} States reached through E-CLOSE(k)
 */
function close(k, nodeClosure, nClosed, nfa) {
    // If E-CLOSE(k) not yet calculated
    if (nodeClosure[k].length == 0) {
        nClosed.add(k);
        // If state k has epsilon transitions
        if (nfa[k].length != 0 && EPSILON in nfa[k]) {
            // For each state immediately reachable via epsilon transitions from state k
            for (var q of nfa[k][EPSILON]) {
                // If state q not in E-CLOSE(n)
                if (!nClosed.has(q)) { 
                    // Add each state from E-CLOSE(q) to E-CLOSE(n)
                    for (var p of close(q, nodeClosure, nClosed, nfa)) { 
                        nClosed.add(p);
                    }
                }
            }
        }
    } else { // if E-CLOSE(k) already calculated
        // Add each state from E-CLOSE(k) to E-CLOSE(n)
        for (var q of nodeClosure[k]) { 
            nClosed.add(q);
        }
    }
    return nClosed.values();
}

/**
 * Extract symbols from transition label
 * @param {String} label Label of transition
 * @returns {IterableIterator<String>} Symbols extracted from 'label'
 */
function getSymbols(label) {
    var s = new Set();
    // Alphabet of NFA
    const alphabet = [];
    for (var c of SIGMA) {
        alphabet.push(c);
    }
    alphabet.push(EPSILON);
    // If any character of label is in alphabet, add to 's'
    for (var char of label) {
        if (alphabet.includes(char)) {
            s.add(char);
        }
    }
    return s.values();
}

/**
 * Create state transition table and accept states from
 *  user defined NFA
 * @returns {
 *  Array<{
 *      table: 
 *          Array<{
 *              stateID: Number,
 *              symbol: String,
 *              stateIDs: Array<Number>
 *          }>,
 *      accept: Array<Number>
 *  }>
 * } State transition table and accept states
 */
function transTable() {
    // State transition table
    const nfa = {};
    // Accept states
    const final = [];
    // Alphabet of NFA
    const alphabet = [];
    for (var c of SIGMA) {
        alphabet.push(c);
    }
    alphabet.push(EPSILON);
    // For each state in user NFA
    for (var n of nodes) {
        // Add state to state transition table
        nfa[n.id] = {};
        // Initialise transitions
        for (var s of alphabet) {
            nfa[n.id][s] = [];
        }
        // Add to array of accept states if applicable
        if (n.accept) {
            final.push(n.id);
        }
    }
    // For each transition in user NFA
    for (var e of edges) {
        // Get valid characters from label
        var symbols = getSymbols(e.label);
        // Add each transition to state transition table
        for (var s of symbols) {
            nfa[e.fromNode.id][s].push(e.toNode.id);
        }
    }

    return {
        table : nfa,
        accept : final
    }
}

/**
 * Construct a circle given three points
 * @param {Number} x1 x-coordinate of point 1
 * @param {Number} y1 y-coordinate of point 1
 * @param {Number} x2 x-coordinate of point 2
 * @param {Number} y2 y-coordinate of point 2
 * @param {Number} x3 x-coordinate of point 3
 * @param {Number} y3 y-coordinate of point 3
 * @returns {
 *  Array<{
 *      x: Number,
 *      y: Number,
 *      radius: Number
 *  }>
 * } Coordinates for centre of circle and the circle's radius
 */
function circleFromPoints(x1, y1, x2, y2, x3, y3) {
    // Find circle equation from three points (above)
    var a = x1*(y2-y3)-y1*(x2-x3)+x2*y3-x3*y2;
    var b = (x1**2+y1**2)*(y3-y2)+(x2**2+y2**2)*(y1-y3)+(x3**2+y3**2)*(y2-y1);
    var c = (x1**2+y1**2)*(x2-x3)+(x2**2+y2**2)*(x3-x1)+(x3**2+y3**2)*(x1-x2);

    var x = -b/(2*a); // x centre
    var y = -c/(2*a); // y centre

    return {
        'x' : x,
        'y' : y,
        'radius' : Math.hypot(x-x1, y-y1)
    }
}

/**
 * Return index from state/edge array given a state/edge ID respectively
 * @param {Number} id State/Edge ID
 * @param {Number} arr State/Edge array
 * @returns {Number} Index of state/edge
 */
function getFromId(id, arr) {
    for (var i=0; i<arr.length; i++) {
        if (arr[i].id == id) {
            return i;
        }
    }
}

/**
 * Return index of edge at the cursor's position
 *  (or -1 if cursor isn't hovering over any edge)
 * @param {Number} x x-position of cursor
 * @param {Number} y y-position of cursor
 * @returns {Number} Index of edge, or -1 if no edge present
 */
function edgeUnderMouse(xm, ym) {
    // For each edge
    for (var i=edges.length-1; i >=0; i--) {
        var edge = edges[i];
        if (edge.id != startTid) { // ignore start edge
            if (edge.fromNode == edge.toNode) { // self-loop
                var dx = edge.x-xm;
                var dy = edge.y-ym;
                // If cursor contained within area of self-loop
                if (dx*dx+dy*dy < (edge.radius+SELECTAREA)*(edge.radius+SELECTAREA)) {
                    return i;
                }
            } else { // normal transition
                var xf = edge.fromNode.x;
                var yf = edge.fromNode.y;
                var xt = edge.toNode.x;
                var yt = edge.toNode.y;
                var dx = xt - xf;
                var dy = yt - yf;
                var len = Math.sqrt(dx*dx+dy*dy);
                if (edge.curved) { // curved transition
                    var ang = 1.5*Math.PI-edge.angle;
                    var cosShift = 2*SELECTAREA*Math.cos(ang);
                    var sinShift = 2*SELECTAREA*Math.sin(ang);
                    // Shifts selection window accordingly
                    var perc = (dx*(xm-xf+cosShift)+dy*(ym-yf-sinShift))/(len*len);
                    var dist = (dx*(ym-yf-sinShift)-dy*(xm-xf+cosShift))/len;
                } else { // straight transition
                    // Regular selection window for transition
                    var perc = (dx*(xm-xf)+dy*(ym-yf))/(len*len); // how far along transition the mouse is (parallel proportion)
                    var dist = (dx*(ym-yf)-dy*(xm-xf))/len; // perpendicular distance mouse is away from line equation of transition
                }
                // If cursor within selection area
                if (perc > 0 && perc < 1 && Math.abs(dist) < SELECTAREA) {
                    return i;
                }
            }
        }
    }
    // If no edge is under cursor
    return -1;
}

/**
 * Return index of node at the cursor's position
 *  (or -1 if cursor isn't hovering over any node)
 * @param {Number} x x-position of cursor
 * @param {Number} y y-position of cursor
 * @returns {Number} Index of node, or -1 if no node present
 */
function nodeUnderMouse(x, y) {
    for (var i=nodes.length-1; i >= 0; i--) {
        var node = nodes[i];
        var dx = node.x-x;
        var dy = node.y-y;
        // Use Pythagoras' Theorem to check if mouse is within node's area
        if (dx*dx+dy*dy < RADIUS*RADIUS) {
            return i;
        }
    }
    return -1 // indicates no node under mouse
}

/**
 * Calculate offset of mouse coordinates from the position of the canvas
 * @param {MouseEvent} event Event when mouse interacts with HTML document
 * @returns {
 *  Array<{
 *      x: Number,
 *      y: Number
 *  }>
 * } Coordinates of mouse with respect to canvas position
 */
function coordinates(event) {
    var dimensions = canvas.getBoundingClientRect();
    // Account for canvas offset by subtracting its top most- and left most-position in window
    return {
        x: event.clientX-dimensions.left,
        y: event.clientY-dimensions.top
    }
}

/**
 * Draws edges and nodes to the canvas
 * @param {Boolean} mouseDown Indicates if the mouse is pressed down
 */
function updateCanvas(mouseDown) {
    // Only update canvas if user is dragging state, pressing key, or clicking mouse
    if (state && (state.dragging || mouseDown)) {
        ctx.clearRect(0, 0, canvas.width, canvas.height); // clear canvas

        // Draw edges
        for (var i=0; i<edges.length; i++) {
            if (edges[i].id != startTid) {
                edges[i].draw(ctx);
            }
        }

        // Draw nodes
        for (var j=0; j<nodes.length; j++) {
            nodes[j].draw(ctx);
            // Draw start edge
            if (nodes[j].id == startSid) {
                edges[getFromId(startTid, edges)].draw(ctx);
            }
        }
    }
}

/**
 * Returns true if the state transition table is one for a DFA
 * @param {
 *  Array<{
 *      stateID: Number,
 *      symbol: String,
 *      stateIDs: Array<Number>
 *  }>
 * } table state transition table for an NFA
 * @returns {Boolean} True if transition table is one for a DFA
 */
function dfaTest(table) {
    for (const symbols of Object.values(table)) {
        for (const [symbol, ids] of Object.entries(symbols)) {
            if (symbol != EPSILON) {
                if (ids.length != 1) {
                    return false;
                }
            } else {
                if (ids.length != 0) {
                    return false;
                }
            }
        }
    }
    return true;
}

/**
 * Outputs correct if regex and user's machine accept the same language
 */
function comp() {
    if (startSid != -1) {
        var proceed = true;
        var userNFA = transTable();
        if (document.getElementById('dfa-toggle').checked) {
            proceed = dfaTest(userNFA.table);
        }
        if (proceed) {
            var regNFA = regularExpression.nfa;
            var user = subsetConstruct(userNFA.table, startSid, userNFA.accept, 0);
            var reg = subsetConstruct(regNFA.table, regNFA.start, [regNFA.end], Object.keys(user.dfa).length)
            var equal = isomorphic(user, reg);
            if (equal) {
                answer.innerHTML = "Correct";
                answer.style.color = "green";
            } else {
                answer.innerHTML = "Incorrect";
                answer.style.color = "red";
            }
        } else {
            answer.innerHTML = "Not a DFA";
            answer.style.color = "red";
        }
    } else {
        answer.innerHTML = "No Start State";
        answer.style.color = "red";
    }
}

/**
 * Generate new regular expression and display to user
 */
function gen() {
    regularExpression.generate();
    regex.innerHTML = regularExpression.regex;
    answer.innerHTML = "Draw A Machine";
    answer.style.color = "black";
}

const canvas = document.getElementById('flat-canvas');
const regex = document.getElementById('regex');
const answer = document.getElementById('result');
const ctx = canvas.getContext('2d');
ctx.fillStyle = STATEFILL;
ctx.textAlign = "center";
ctx.font = FONTSIZE + "px Arial";
var fromX = 0;
var fromY = 0;

const regularExpression = new Regex();
regex.innerHTML = regularExpression.regex;
answer.innerHTML = "Draw A Machine";

var state = null;

/**
 * Pressing a key
 */
window.addEventListener("keydown",
    function(event){

        // State/Edge to add label to
        var addLabel = null;

        if (highSid != -1) { // state highlighted
            index = getFromId(highSid, nodes);
            addLabel = nodes[index];
        } else if (highTid != -1) { // edge highlighted
            index = getFromId(highTid, edges);
            addLabel = edges[index];
        }

        // If State/Edge highlighted and character entered
        if (addLabel != null && event.key.length == 1) {
            // Length of current label
            var length = addLabel.label.length;
            // Convert '\e' into EPSILON if read, else add character
            if (length > 0 && length < 7 && addLabel.label[length-1] == '\\' && event.key == 'e') {
                addLabel.label = addLabel.label.slice(0,-1) + EPSILON;
            } else if (length < 6) {
                addLabel.label += event.key;
            }
        } else if (event.key == "Backspace") { // delete last character
            addLabel.label = addLabel.label.slice(0,-1);
        } else if (event.key == "Delete") { // delete highlighted State/Edge
            if (highSid != -1) { // state selected
                // Create new edge set
                var new_edges = [];
                // Only include edges not adjacent to selected node
                for (var i=0; i<edges.length; i++) {
                    if (edges[i].fromNode == nodes[index] || edges[i].toNode == nodes[index]) {
                        if (edges[i].id == startTid) {
                            startTid = -1;
                        }
                    } else {
                        new_edges.push(edges[i]);
                    }
                }
                edges = new_edges;
                if (nodes[index].id == startSid) {
                    startSid = -1;
                }
                // Remove state
                nodes.splice(index,1);
                highSid = -1;
            } else if (highTid != -1) { // edge selected
                if (edges[index].id == startTid) {
                    startTid = -1;
                }
                for (var i=0; i<edges.length; i++) {
                    // If selected edge curved, remove curved property of other curved edge
                    if (edges[i].fromNode == edges[index].toNode && edges[i].toNode == edges[index].fromNode) {
                        edges[i].curved = false;
                        break;
                    }
                }
                edges.splice(index,1);
                highTid = -1;
            }
        } else if (event.key == "Escape") {
            highSid = -1;
            highTid = -1;
        }

        updateCanvas(true);
    }
);

/**
 * Double-clicking on the canvas
 */
canvas.addEventListener("dblclick",
    function(event) {

        // Get mouse coordinates
        var coords = coordinates(event);
        // x- and y-coordinates of mouse
        var x = coords.x;
        var y = coords.y;
        // Get State/Edge under mouse
        var stateIndex = nodeUnderMouse(x, y);
        var edgeIndex = edgeUnderMouse(x, y);

        if (stateIndex != -1) { // node selected
            // Toggle accept state
            nodes[stateIndex].accept = !nodes[stateIndex].accept;
        } else if (edgeIndex == -1) { // empty space on canvas selected
            if (event.shiftKey) { // shift held
                if (highSid != -1) { // state highlighted
                    // Create new node
                    var n = new Node(sid, x, y);
                    state = n;
                    nodes.push(n);
                    // Add edge between highlighted state and new state
                    var e = new Edge(tid, nodes[getFromId(highSid, nodes)], n);
                    sid++;
                    edges.push(e);
                    tid++;
                }
            } else if (event.ctrlKey) { // ctrl held
                // Create new node
                var n = new Node(sid, x, y);
                state = n;
                nodes.push(n);
                highSid = sid;
                highTid = -1;
                startSid = sid;
                sid++;
                // Create start edge if it doesn't exist, and point it to new node
                if (startTid == -1) {
                    var e = new Edge(tid, null, n);
                    edges.push(e);
                    startTid = tid;
                    tid++;
                } else {
                    // Set existing start edge to point at this node
                    for (var i=0; i<edges.length; i++) {
                        if (edges[i].id == startTid) {
                            edges[i].toNode = nodes[getFromId(highSid, nodes)];
                            break;
                        }
                    }
                }
            } else { // no extra key held
                // Create new node and highlight it
                var n = new Node(sid, x, y);
                state = n;
                nodes.push(n);
                highSid = sid;
                highTid = -1;
                sid++;
            }
        }
        updateCanvas(true);
    }
);

/**
 * Mouse clicked
 */
canvas.addEventListener("mousedown",
    function(event) {

        // Get mouse coordinates
        var coords = coordinates(event);
        // x- and y-coordinates of mouse
        var x = coords.x;
        var y = coords.y;
        // Get State/Edge under mouse
        var stateIndex = nodeUnderMouse(x, y);
        var edgeIndex = edgeUnderMouse(x, y);

        if (stateIndex != -1) { // state selected
            if (event.shiftKey) { // shift held
                if (highSid != -1) { // a state is currently highlighted
                    var from = nodes[getFromId(highSid, nodes)];
                    var create = true;
                    var curve = false;
                    for (var i=0; i<edges.length; i++) {
                        // If edge already exists between nodes, do not create a new one
                        if (create && edges[i].fromNode == from && edges[i].toNode == nodes[stateIndex]) {
                            create = false;
                        }
                        // If reversed edge exists, curve both
                        if (!curve && edges[i].fromNode == nodes[stateIndex] && edges[i].toNode == from) {
                            curve = true;
                            edges[i].curved = true;
                        }
                    }
                    if (create) {
                        // Add new edge between nodes
                        var e = new Edge(tid, nodes[getFromId(highSid, nodes)], nodes[stateIndex]);
                        edges.push(e);
                        tid++;
                        e.curved = curve;
                    }
                }
            } else if (event.ctrlKey) { // ctrl held
                state = nodes[stateIndex];
                highSid = state.id;
                highTid = -1;
                startSid = highSid; // set highlighted state as start state

                if (startTid == -1) {
                    // Create start edge if it doesn't exist, and point it to node
                    var e = new Edge(tid, null, state);
                    edges.push(e);
                    startTid = tid;
                    tid++;
                } else {
                    // Set start edge to point at this node
                    for (var i=0; i<edges.length; i++) {
                        if (edges[i].id == startTid) {
                            edges[i].toNode = nodes[getFromId(highSid, nodes)];
                            break;
                        }
                    }
                }
            } else { // no extra key held
                // User drags node
                state = nodes[stateIndex];
                state.dragging = true;
                highSid = state.id;
                highTid = -1;
                canvas.style.cursor = "move";
            }
        } else if (edgeIndex != -1) { // edge selected
            // Highlight edge
            var edge = edges[edgeIndex];
            highTid = edge.id;
            highSid = -1;
        } 
        updateCanvas(true);
    }
);

/**
 * Mouse moving over canvas
 */
canvas.addEventListener("mousemove",
    function(event) {
        var coords = coordinates(event); // get mouse coordinates
        var x = coords.x;
        var y = coords.y;
        var stateId = nodeUnderMouse(x,y);

        // Change look of mouse if hovering over state
        if (stateId != -1) {
            canvas.style.cursor = "move";
        } else {
            canvas.style.cursor = "auto";
        }

        // Calculate change in mouse position
        var dx = x-fromX;
        var dy = y-fromY;
        fromX = x;
        fromY = y;

        if (state && state.dragging) { // there exists a state being dragged
            state.x += dx;
            state.y += dy;
            updateCanvas(false); // only update if dragging node
        }
    }
);

/**
 * Mouse click released
 */
canvas.addEventListener("mouseup",
    function(){
        // If canvas nonempty and dragging state
        if (state && state.dragging) {
            // Stop dragging state
            state.dragging = false;
            // Return cursor style to pointer
            canvas.style.cursor = "auto";
        }
    }
);
