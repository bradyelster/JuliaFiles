using PkgTemplates

t = Template(;
    user = "bradyelster",
    dir = "/home/belster/TearingModes/",
    julia = v"1.11",
    authors = ["Brady Elster"],
    plugins = [
        License(; name = "MIT"),
        Logo(),
        TagBot(),
        GitHubActions(; destination = "CI.yml", coverage = true),
        Codecov(),
        Secret("MySecret"),
        # RegisterAction(),
        Formatter(; style = "sciml"),
        Git(; branch = "main", jl = true, manifest = true, ssh = true),
        Documenter{GitHubActions}(),
        Develop(),
    ]
)

generate("TearingModes", t)
