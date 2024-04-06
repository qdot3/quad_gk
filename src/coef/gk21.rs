use crate::QuadGKCoef;

/// Shortcut to the nodes and weights of 21st order Gauss Kronrod quadrature rule.
pub struct GK21;

impl QuadGKCoef<21> for GK21 {
    fn quad_gk_coef(&self) -> [(f64, f64, Option<f64>); 21] {
        [
            (
                -0.9956571630258080807355272806890029,
                0.011694638867371874278064396062192,
                None,
            ),
            (
                -0.9739065285171717200779640120844521,
                0.0325581623079647274788189724593898,
                Some(0.066671344308688137593568809893332),
            ),
            (
                -0.9301574913557082260012071800595084,
                0.0547558965743519960313813002445802,
                None,
            ),
            (
                -0.8650633666889845107320966884234931,
                0.07503967481091995276704314091619001,
                Some(0.149451349150580593145776339657697),
            ),
            (
                -0.7808177265864168970637175783450424,
                0.09312545458369760553506546508336634,
                None,
            ),
            (
                -0.6794095682990244062343273651148736,
                0.109387158802297641899210590325805,
                Some(0.2190863625159820439955349342281632),
            ),
            (
                -0.5627571346686046833390000992726941,
                0.1234919762620658510779581098310742,
                None,
            ),
            (
                -0.4333953941292471907992659431657842,
                0.1347092173114733259280540017717068,
                Some(0.2692667193099963550912269215694694),
            ),
            (
                -0.2943928627014601981311266031038656,
                0.1427759385770600807970942731387171,
                None,
            ),
            (
                -0.14887433898163121088482600112972,
                0.1477391049013384913748415159720681,
                Some(0.2955242247147528701738929946513383),
            ),
            (0.0, 0.1494455540029169056649364683898212, None),
            (
                0.14887433898163121088482600112972,
                0.147739104901338491374841515972068,
                Some(0.2955242247147528701738929946513383),
            ),
            (
                0.2943928627014601981311266031038656,
                0.142775938577060080797094273138717,
                None,
            ),
            (
                0.4333953941292471907992659431657842,
                0.134709217311473325928054001771707,
                Some(0.2692667193099963550912269215694694),
            ),
            (
                0.5627571346686046833390000992726941,
                0.123491976262065851077958109831074,
                None,
            ),
            (
                0.6794095682990244062343273651148736,
                0.109387158802297641899210590325805,
                Some(0.2190863625159820439955349342281632),
            ),
            (
                0.7808177265864168970637175783450424,
                0.093125454583697605535065465083366,
                None,
            ),
            (
                0.8650633666889845107320966884234931,
                0.07503967481091995276704314091619,
                Some(0.1494513491505805931457763396576973),
            ),
            (
                0.9301574913557082260012071800595084,
                0.05475589657435199603138130024458,
                None,
            ),
            (
                0.9739065285171717200779640120844521,
                0.0325581623079647274788189724593898,
                Some(0.0666713443086881375935688098933318),
            ),
            (
                0.9956571630258080807355272806890029,
                0.011694638867371874278064396062192,
                None,
            ),
        ]
    }
}
