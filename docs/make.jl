using Documenter, Kinetica, KineticaKPM

makedocs(
    sitename = "Kinetica Documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold = 307200
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Tutorials" => [
            "Arbitrary Simulation Conditions" => "tutorials/arbitrary-conditions.md",
            "ODE Solution" => "tutorials/ode-solution.md",
            "Kinetic Calculators" => "tutorials/kinetic-calculators.md",
            "Iterative CRN Exploration" => "tutorials/iterative-exploration.md",
            "Results Analysis" => "tutorials/results-analysis.md",
            "Filtering CRNs" => "tutorials/filtering-crns.md",
            "Saving & Loading" => "tutorials/saving-loading.md",
            "Logging" => "tutorials/logging.md"
        ],
        "Developing with Kinetica" => [
            "CRN Representation" => "development/crn-representation.md",
            "Condition Profiles" => "development/condition-profiles.md",
            "Calculator Interface" => "development/calculator-interface.md"
        ],
        "API" => [
            "Kinetica.jl" => [
                "Exploration" => "api/kinetica/exploration.md",
                "Solving" => "api/kinetica/solving.md",
                "Conditions" => "api/kinetica/conditions.md",
                "Analysis" => "api/kinetica/analysis.md",
                "Open Babel" => "api/kinetica/openbabel.md",
                "RDKit" => "api/kinetica/rdkit.md",
                "Utilities" => "api/kinetica/utilities.md"
            ],
            "KineticaKPM.jl" => "api/kineticakpm.md"
        ]
    ]
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/Kinetica-jl/Kinetica.jl.git",
        push_preview = true
    )
end