#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PARAM_NUM 40

typedef struct {
    double x;
    double y;
} point;

typedef struct {
    double top;
    double bottom;
} interval;

int n = 10000;
int interval_num = 10000; // должен быть четным
interval *v;
double a[PARAM_NUM];

void init_sofa(int interval_num, interval x[]) {
    for (int i = interval_num + 1; i--;) {
        x[i].top = 1;
        x[i].bottom = 0;
    }
}

point get_lower_corner(const double alpha, const double params[]) {
    point p;

    p.x = 0;// params[0] * cos(2 * alpha);
    p.y = 0;// params[0] * sin(2 * alpha);

    for (int j = 0; j < PARAM_NUM / 2; j++) {
        p.x += params[2 * j] * cos(alpha * 2 * (2 * j + 1));
        p.y += params[2 * j + 1] * sin(alpha * 2 * (2 * j + 1));
    }

    return p;
}

point get_upper_corner(const double alpha, const point lower_corner) {
    point p;
    p.x = lower_corner.x + sqrt(2) * cos(M_PI / 4 + alpha);
    p.y = lower_corner.y + sqrt(2) * sin(M_PI / 4 + alpha);

    return p;
}

// Интегрирование методом Симпсона
double integrate(const double sofa_length, const int interval_num, const interval *x) {
    double result = 0.0;

    if (x[0].top > x[0].bottom) {
        result += (x[0].top - x[0].bottom);
    }

    int j = 1;
    while (j < interval_num) {
        if (x[j].top > x[j].bottom) {
            result += 4 * (x[j].top - x[j].bottom);
        }
        j++;

        if (x[j].top > x[j].bottom) {
            result += 2 * (x[j].top - x[j].bottom);
        }
        j++;
    }

    // Выше в цикле прибавили лишнего, вычитаем вместо прибавления
    if (x[interval_num].top >= x[interval_num].bottom) {
        result -= (x[interval_num].top - x[interval_num].bottom);
    }

    return result / interval_num * sofa_length / 3.0;
}

double get_area(const int interval_num, const int positions_num, const double *params, interval *x) {
    const double sofa_length = 2 * get_lower_corner(0, params).x + 2;

    init_sofa(interval_num, x);

    for (int i = 1; i < positions_num; i++) {
        int j, k;
        double y, shift_at_zero_index, boundary_change_for_interval;

        double alpha = M_PI * (double) i / positions_num * 0.5;// + PI/positions_num*0.25;
        double tg = tan(alpha);
        double ctg = -1 / tg;

        point lower_corner = get_lower_corner(alpha, params);
        point upper_corner = get_upper_corner(alpha, lower_corner);

        // Top curve
        k = (int) ((upper_corner.x / sofa_length + 0.5) * interval_num);
        boundary_change_for_interval = tg * sofa_length / interval_num;
        shift_at_zero_index = upper_corner.y - tg * (upper_corner.x + sofa_length * 0.5);
        y = shift_at_zero_index;
        for (j = 0; j < k; j++) {
            if (x[j].top > y)
                x[j].top = y;
            y += boundary_change_for_interval;
        }

        k = (int) ((upper_corner.x / sofa_length + 0.5) * interval_num);
        boundary_change_for_interval = ctg * sofa_length / interval_num;
        shift_at_zero_index = upper_corner.y - ctg * (upper_corner.x + sofa_length * 0.5);
        y = shift_at_zero_index + boundary_change_for_interval * (k + 1);
        for (j = k + 1; j < interval_num + 1; j++) {
            if (x[j].top > y)
                x[j].top = y;
            y += boundary_change_for_interval;
        }

        // Bottom curve
        k = (int) ((lower_corner.x / sofa_length + 0.5) * interval_num);
        boundary_change_for_interval = tg * sofa_length / interval_num;
        shift_at_zero_index = lower_corner.y - tg * (lower_corner.x + sofa_length * 0.5);
        y = shift_at_zero_index;
        for (j = 0; j < k; j++) {
            if (x[j].bottom < y)
                x[j].bottom = y;
            y += boundary_change_for_interval;
        }
        boundary_change_for_interval = ctg * sofa_length / interval_num;
        shift_at_zero_index = lower_corner.y - ctg * (lower_corner.x + sofa_length * 0.5);
        y = shift_at_zero_index + boundary_change_for_interval * (k + 1);
        for (j = k + 1; j < interval_num + 1; j++) {
            if (x[j].bottom < y)
                x[j].bottom = y;
            y += boundary_change_for_interval;
        }
    }

    return integrate(sofa_length, interval_num, x);
}

void maximize(const int interval_num, const int positions_num, double params[], interval x[]) {
    double grad[PARAM_NUM], new_params[PARAM_NUM], q = 0.5;
    int i;

    // С чего-то начинаем
    double f = get_area(interval_num, positions_num, params, x);
    double f_prev = f - 0.1;
    double f_prev2 = f - 0.2;

    while (1 /*positions_num*positions_num*fabs(f-f_prev) > 0.01*(1-q)*/) {

        // Вычисляем градиент
        for (i = 0; i < PARAM_NUM; i++) {
            double old_param = params[i];

            params[i] += 1.0e-9;
            double changed_f = get_area(interval_num, positions_num, params, x);
            params[i] = old_param;

            grad[i] = (changed_f - f) * 1.0e9;
//			printf("i=%d; int = %f; grad = %f; \positions_num", i, changed_f, grad[i]);
        }

        // Пытаемся понять, на какой шаг можно пройти вперед.

        // Сначала отступаем немного.
        double lambda = 1.0e-5;
        for (i = 0; i < PARAM_NUM; i++) {
            new_params[i] = params[i] + lambda * grad[i];
        }
        double f_at_lambda = get_area(interval_num, positions_num, new_params, x);

        printf("               int = %5.10f; lambda=%e\n", f_at_lambda, lambda);

        double multiplier = (f_at_lambda > f) ? 10 : 0.1;
        while (lambda < 100 && lambda > 1.0e-10) {
            // Затем отступаем больше и сравниваем, что получилось
            double lambda2 = lambda * multiplier;
            for (i = 0; i < PARAM_NUM; i++) {
                new_params[i] = params[i] + lambda2 * grad[i];
            }
            double f_at_lambda2 = get_area(interval_num, positions_num, new_params, x);

            printf("               int = %5.10f; lambda=%e\n", f_at_lambda2, lambda2);

            // Увеличение шага ничего не дает, ну и ладно.
            if (f_at_lambda2 < f_at_lambda && f_at_lambda > f) {
                // Если шаг уменьшаем, убеждаемся, что действительно нашли точку максимума лучше, чем было
                break;
            }

            // Увеличение помогло, пробуем еще.
            lambda = lambda2;
            f_at_lambda = f_at_lambda2;
        }

        for (i = 0; i < PARAM_NUM; i++) {
            // Отступаем на найденный шаг по найденному градиенту.
            params[i] += lambda * grad[i];
            printf("a[%i] = %5.80f;\n", i, params[i]);
        }

        f_prev2 = f_prev;
        f_prev = f;
        f = f_at_lambda;
        q = (f - f_prev) / (f_prev - f_prev2);
        printf("-------------->>> int = %5.10f; q=%5.10f; lambda=%e, positions_num=%d\n", f, q, lambda, positions_num);

        if (lambda <= 1.0e-10) {
            break;
        }

        if (q < 0) {
            printf("ERROR\n");
            exit(0);
        }
        if (q < 0.5) {
            q = 0.5;
        }
    }

}

int main() {

/*	a[0] = 0.602249;
	a[1] = 0.004820;
	a[2] = 0.664509;// 0.63661977236758;
	a[3] = 0.006013;
	a[4] = 0.001659;
	a[5] = 0.000422;
	a[6] = 0.000249;
	a[7] = 0.000058;
	a[8] = 0.001329;
	a[9] = 0.000630;
*/
    a[0] = 0.60265324154425170544158163465908728539943695068359375000000000000000000000000000;
    a[1] = 0.66507963190437369149776714039035141468048095703125000000000000000000000000000000;
    a[2] = 0.00453773100691790172434014749569541891105473041534423828125000000000000000000000;
    a[3] = 0.00588952547592623033845260493990281247533857822418212890625000000000000000000000;
    a[4] = 0.00153583989475173429001264580051611119415611028671264648437500000000000000000000;
    a[5] = 0.00175907100299416228189608446541569719556719064712524414062500000000000000000000;
    a[6] = 0.00071754172280549105827213152153376540809404104948043823242187500000000000000000;
    a[7] = 0.00045681566437125496011814607122403231187490746378898620605468750000000000000000;
    a[8] = 0.00041501851707649519036824603546165235457010567188262939453125000000000000000000;
    a[9] = 0.00029094140850572568430046360710150565864751115441322326660156250000000000000000;
    a[10] = 0.00031011609052774332523361167091024981345981359481811523437500000000000000000000;
    a[11] = 0.00012956069822851217139027257818639782271930016577243804931640625000000000000000;
    a[12] = 0.00019478031459025047173252986887348470190772786736488342285156250000000000000000;
    a[13] = 0.00009410431256871465042438545944847305690927896648645401000976562500000000000000;
    a[14] = 0.00016046766343685628093511386094149884229409508407115936279296875000000000000000;
    a[15] = 0.00005792883490757016811886770391026857396354898810386657714843750000000000000000;
    a[16] = 0.00010772002913539858049568609388302320439834147691726684570312500000000000000000;
    a[17] = 0.00003778814825940114609931219646909994480665773153305053710937500000000000000000;
    a[18] = 0.00008645587684857233248987562479470625476096756756305694580078125000000000000000;
    a[19] = 0.00003519308996982051376345781990018224405503133311867713928222656250000000000000;
    a[20] = 0.00006318477054875420799960605844347583115450106561183929443359375000000000000000;
    a[21] = 0.00001984241807511111165235155595620852864158223383128643035888671875000000000000;
    a[22] = 0.00004431844865500079696032931231286511319922283291816711425781250000000000000000;
    a[23] = 0.00001631946427648855906856507902435993173639872111380100250244140625000000000000;
    a[24] = 0.00003503726635226462227286348127819337605615146458148956298828125000000000000000;
    a[25] = 0.00001035462496703942469617864413633512299384165089577436447143554687500000000000;
    a[26] = 0.00001231895210328293073498689275124817754658579360693693161010742187500000000000;
    a[27] = 0.00001358541624575160784520709567058105449177674017846584320068359375000000000000;
    a[28] = 0.00001452017027942621489545928920028572406408784445375204086303710937500000000000;
    a[29] = 0.00001053854992125293686634421308490061619522748515009880065917968750000000000000;


    for (n = 1000; n < 1000000; n *= 10) {
        interval_num = n;
        printf("----------------->>> n = %d, interval_num = %d\n", n, interval_num);
        v = (interval *) malloc((interval_num + 1) * sizeof(interval));
        maximize(interval_num, n, a, v);
        free(v);
    }
//	printf("int = %f", get_area());
    return 0;
}
