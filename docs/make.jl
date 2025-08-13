using Documenter, ChargeFlipPhaser

makedocs(
    sitename = "ChargeFlipPhaser",
    modules = [ChargeFlipPhaser],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api/public.md"
    ]
)

deploydocs(
    repo = "github.com/pavel-kalouguine/ChargeFlipPhaser.jl.git"
)