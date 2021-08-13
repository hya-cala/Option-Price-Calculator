function fbs() {
    var cop = document.getElementById("cop").selectedIndex;
    var S = parseFloat(document.getElementById("S").value);
    var X = parseFloat(document.getElementById("L").value);
    var r = parseFloat(document.getElementById("r").value);
    var v = parseFloat(document.getElementById("vol").value);
    var T = parseFloat(document.getElementById("t").value);

    //ŷʽ��Ȩ��
    if (cop < 2) {
        var d1, d2;

        d1 = (Math.log(S / X) + (r + v * v / 2.0) * T) / (v * Math.sqrt(T));
        d2 = d1 - v * Math.sqrt(T);

        if (cop == 0)
            var result = S * CND(d1) - X * Math.exp(-r * T) * CND(d2);
        else
            var result = X * Math.exp(-r * T) * CND(-d2) - S * CND(-d1);
        document.getElementById("output").innerHTML = "Option Price: " + result + "</br>";
    }
    else {
        if (cop == 2) {
            var N = 2000;
            var t = T / N;
            var u = Math.exp(v * Math.sqrt(t));
            var d = 1 / u;
            var a = Math.exp(r * t);
            var p = (a - d) / (u - d);

            var SM = new Array();
            for (var k = 0; k <= (N + 1); k++) {
                SM[k] = new Array();
                for (var j = 0; j <= (N + 1); j++) {
                    SM[k][j] = 0;
                }
            }

            for (i = 0; i <= N; i++) {
                for (j = 0; j <= i; j++) {
                    SM[i + 1][j + 1] = S * (u ** j) * (d ** (i - j));
                }
            }

            var VM = new Array()
            for (i = 0; i <= N; i++) {
                VM[i + 1] = Math.max(-SM[N + 1][i + 1] + X, 0);
            }

            var Vn = new Array()
            for (j = 1; j <= N; j++) {
                for (i = 1; i <= (VM.length - 1); i++) {
                    St = Math.max(X - SM[N + 1 - j][i], 0);
                    no = Math.exp(-r * t) * ((1 - p) * VM[i] + p * VM[i + 1]);
                    if (St >= no) {
                        Vn[i] = St;
                    }
                    else {
                        Vn[i] = no;
                    }
                }
                VM = Vn;
                Vn = new Array();
            }
            document.getElementById("output").innerHTML = "Option Price: " + VM[1] + "</br>";
        }
        else {
            var N = 2000;
            var t = T / N;
            var u = Math.exp(v * Math.sqrt(t));
            var d = 1 / u;
            var a = Math.exp(r * t);
            var p = (a - d) / (u - d);

            var SM = new Array();
            for (var k = 0; k <= (N + 1); k++) {
                SM[k] = new Array();
                for (var j = 0; j <= (N + 1); j++) {
                    SM[k][j] = 0;
                }
            }

            for (i = 0; i <= N; i++) {
                for (j = 0; j <= i; j++) {
                    SM[i + 1][j + 1] = S * (u ** j) * (d ** (i - j));
                }
            }

            var VM = new Array()
            for (i = 0; i <= N; i++) {
                VM[i + 1] = Math.max(SM[N + 1][i + 1] - X, 0);
            }

            var Vn = new Array()
            for (j = 1; j <= N; j++) {
                for (i = 1; i <= (VM.length - 1); i++) {
                    St = Math.max(SM[N + 1 - j][i] - X, 0);
                    no = Math.exp(-r * t) * ((1 - p) * VM[i] + p * VM[i + 1]);
                    if (St >= no) {
                        Vn[i] = St;
                    }
                    else {
                        Vn[i] = no;
                    }
                }
                VM = Vn;
                Vn = new Array();
            }
            document.getElementById("output").innerHTML = "Option Price: " + VM[1] + "</br>";
        }
    }

}

/* The cummulative Normal distribution function: */

function CND(x) {
    var a1, a2, a3, a4, a5, k;

    a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937, a4 = -1.821255978, a5 = 1.330274429;

    if (x < 0.0)
        return 1 - CND(-x);
    else
        k = 1.0 / (1.0 + 0.2316419 * x);
    return 1.0 - Math.exp(-x * x / 2.0) / Math.sqrt(2 * Math.PI) * k
        * (a1 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + k * 1.330274429))));
}
