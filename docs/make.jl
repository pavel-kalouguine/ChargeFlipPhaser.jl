using Documenter, ChargeFlipPhaser

makedocs(
    sitename = "ChargeFlipPhaser",
    modules = [ChargeFlipPhaser],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api/public.md",
        "Internal API" => "api/private.md"
    ]
)

deploydocs(
    repo = "github.com/pavel-kalouguine/ChargeFlipPhaser.jl.git"
)