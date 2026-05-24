{{fullname | escape | underline}}


{% if exists("extra/" ~ module ~ "." ~ objname ~ ".rst") %}
.. include:: ../extra/{{module}}.{{objname}}.rst
{% endif %}

.. currentmodule:: {{module}}

.. autoclass:: {{objname}}

{% block methods %}
    {% if methods %}
    .. rubric:: {{_('Methods')}}

    .. autosummary::
        :nosignatures:
        {% for item in methods %}
            {%- if not item.startswith('_') or item.endswith('__') %}
        {{name}}.{{item}}
            {%- endif -%}
        {%- endfor %}
    {% endif %}
{% endblock %}

{% block attributes %}
    {% if attributes %}
    .. rubric:: {{_('Attributes')}}

    .. autosummary::
        {% for item in attributes %}
        {{name}}.{{item}}
        {%- endfor %}
    {% endif %}
{% endblock %}

{% for item in methods %}
    .. automethod:: {{name}}.{{item}}

        {% if exists("extra/" ~ module ~ "." ~ name ~ "." ~ item ~ ".rst") %}
        .. include:: ../extra/{{module}}.{{name}}.{{item}}.rst
        {% endif %}
{%- endfor %}
