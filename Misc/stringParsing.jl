function transform(ch::Char)
    ch == '-' ? ch + 50 :
    isspace(ch) || isdigit(ch) ? "" :
    isuppercase(ch) ? "-" * lowercase(ch) :
    'α' <= ch <= 'ω' ? "?" :
    string(ch)
end

function clean(s::String)
    join(transform(ch) for ch in s)
end