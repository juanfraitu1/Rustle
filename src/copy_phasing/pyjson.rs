//! Python's text forms where a product depends on them: `json.dumps` defaults (the `--dump` lines), `repr(float)`, `'%.3f'` and
//! the `repr` of a str-keyed dict (the progress line).

/// A JSON value in Python's insertion order.
#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(String),
    List(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    pub fn str(s: impl Into<String>) -> Json {
        Json::Str(s.into())
    }

    /// `json.dumps(value)`: separators ", " and ": ", ensure_ascii.
    pub fn dumps(&self) -> String {
        let mut out = String::new();
        self.write(&mut out);
        out
    }

    fn write(&self, out: &mut String) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Int(n) => out.push_str(&n.to_string()),
            Json::Float(x) => out.push_str(&float_repr(*x)),
            Json::Str(t) => write_str(t, out),
            Json::List(xs) => {
                out.push('[');
                for (k, x) in xs.iter().enumerate() {
                    if k > 0 {
                        out.push_str(", ");
                    }
                    x.write(out);
                }
                out.push(']');
            }
            Json::Obj(kv) => {
                out.push('{');
                for (k, (key, x)) in kv.iter().enumerate() {
                    if k > 0 {
                        out.push_str(", ");
                    }
                    write_str(key, out);
                    out.push_str(": ");
                    x.write(out);
                }
                out.push('}');
            }
        }
    }
}

/// A JSON string with `ensure_ascii` escapes.
fn write_str(t: &str, out: &mut String) {
    out.push('"');
    for ch in t.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 || (c as u32) > 0x7f => {
                let mut buf = [0u16; 2];
                for u in c.encode_utf16(&mut buf) {
                    out.push_str(&format!("\\u{:04x}", u));
                }
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

/// Python `repr(float)` for a finite value.
pub fn float_repr(x: f64) -> String {
    let e = format!("{x:e}"); // shortest round-trip digits: "8.181818181818182e-1"
    let (sign, e) = match e.strip_prefix('-') {
        Some(r) => ("-", r.to_string()),
        None => ("", e),
    };
    let (mant, exp) = e.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    let n = digits.len() as i32;
    let body = if !(-4..16).contains(&exp) {
        let m = if n > 1 { format!("{}.{}", &digits[..1], &digits[1..]) } else { digits.clone() };
        format!("{m}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    } else if exp >= 0 {
        if n <= exp + 1 {
            format!("{digits}{}.0", "0".repeat((exp + 1 - n) as usize))
        } else {
            format!("{}.{}", &digits[..(exp + 1) as usize], &digits[(exp + 1) as usize..])
        }
    } else {
        format!("0.{}{digits}", "0".repeat((-exp - 1) as usize))
    };
    format!("{sign}{body}")
}

/// Python `'%.3f' % x`.
pub fn fmt3(x: f64) -> String {
    format!("{x:.3}")
}

/// `repr(dict)` of a {str: int} dict in insertion order: `{'a b': 2, 'c': 1}`.
pub fn dict_repr(xs: &[(&str, usize)]) -> String {
    let inner: Vec<String> = xs.iter().map(|(k, v)| format!("'{k}': {v}")).collect();
    format!("{{{}}}", inner.join(", "))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_repr_matches_python() {
        for (x, want) in [
            (1.0, "1.0"),
            (0.1, "0.1"),
            (1e-5, "1e-05"),
            (123456789012345680.0, "1.2345678901234568e+17"),
            (0.000123, "0.000123"),
            (0.8181818181818182, "0.8181818181818182"),
            (1e16, "1e+16"),
            (1e15, "1000000000000000.0"),
            (0.0001, "0.0001"),
            (2.5e-5, "2.5e-05"),
            (0.42857142857142855, "0.42857142857142855"),
            (0.0, "0.0"),
        ] {
            assert_eq!(float_repr(x), want, "{x}");
        }
    }

    #[test]
    fn fmt3_rounds_like_python() {
        for (x, want) in [(0.0625, "0.062"), (0.0005, "0.001"), (0.1875, "0.188"), (2.0 / 3.0, "0.667"), (0.3125, "0.312"), (1.0, "1.000")] {
            assert_eq!(fmt3(x), want, "{x}");
        }
    }

    #[test]
    fn dumps_matches_json_dumps() {
        let v = Json::Obj(vec![(
            "a".into(),
            Json::List(vec![Json::Int(1), Json::Null, Json::Bool(true), Json::str("x\"y\\"), Json::str("é"), Json::Float(0.5)]),
        )]);
        assert_eq!(v.dumps(), concat!(r#"{"a": [1, null, true, "x\"y\\", ""#, "\\u00e9", r#"", 0.5]}"#));
        assert_eq!(Json::Obj(vec![]).dumps(), "{}");
        assert_eq!(Json::List(vec![]).dumps(), "[]");
    }

    #[test]
    fn dict_repr_matches_python() {
        assert_eq!(dict_repr(&[("unexplained", 2), ("equals a reference copy", 2)]), "{'unexplained': 2, 'equals a reference copy': 2}");
        assert_eq!(dict_repr(&[]), "{}");
    }
}
