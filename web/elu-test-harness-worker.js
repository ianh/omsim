importScripts('elu-test-harness.js');

var omsim;
var out = "";
function output(str) {
    if (!str)
        return;
    out = str + "\n" + out;
}
var queue = []
OMSIM({
    print: output,
    printErr: output,
    noInitialRun: true
}).then(function (result) {
    omsim = result;
    omsim._init();
    for (var i = 0; i < queue.length; ++i)
        handle_message(queue[i]);
    queue = [];
});
function handle_message(m) {
    if (m.type == 'solution') {
        var arr = new Uint8Array(m.solution);
        omsim.ccall('load_solution', null, ['array', 'number'], [arr, arr.length]);
        for (var b = m.start; b < m.end; ++b) {
            for (var a = -128; a < 128; ++a) {
                if (a - b < -128 || a - b > 127)
                    continue;
                out = "";
                var r = -999;
                try {
                    var r = omsim._test(a, b, m.fun);
                } catch (e) {
                    output("exception on cycle " + omsim._last_cycle() + ": " + e);
                    output("capacity=" + omsim._last_capacity() + ", used=" + omsim._last_used() + ", removed=" + omsim._last_removed() + ", query=(" + omsim._last_query_u() + ", " + omsim._last_query_v() + ") on most recent atom hash table call");
                } finally {
                    postMessage({type: 'result', a: a, b: b, r: r, output: out });
                    if (r == -999)
                        return;
                }
            }
        }
    }
}
self.addEventListener('message', function (e) {
    if (omsim)
        handle_message(e.data);
    else
        queue.push(e.data);
});
